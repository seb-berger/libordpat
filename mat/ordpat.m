function n = ordpat(x, ord, lag, algorithm)
%ORDPAT Ordinal pattern extraction and encoding.
%   N = ORDPAT(X, ORD, LAG, ALGORITHM) extracts the consecutive ordinal
%   patterns of order ORD from the time series contained in X, using the
%   time lag LAG and the encoding algorithm ALGORITHM.
%
%   N = ORDPAT(X, ORD, LAG) automatically selects a suitable encoding
%   algorithm ('ordpat_c', 'ordpat_mp_c', or 'vectorised') depending on
%   the value of ORD and the availability of MEX library functions.
%
%   N = ORDPAT(X, ORD) uses a default time lag LAG = 1, and automatically
%   selects a suitable encoding algorithm ('overlap_c' or 'overlap_mp_c')
%   depending on the value of ORD.
%
%   ORDPAT('--make') compiles the underlying C functions into MEX files.
%   This has to be called only once. After successful completion, the
%   pre-compiled encoding algorithms can be used. A correctly configured
%   MEX environment and a supported C compiler are required. If in doubt,
%   see the section on MEX in the GNU Octave or MATLAB user manual.
%
%   ORDPAT('--clean') removes all MEX files created by 'ordpat('--make')'.
%
%   ORDPAT always operates along the first dimension of X, that is, each
%   column of the input matrix is treated as an individual time series.
%
%   The ORDPAT function provides various encoding algorithms, which all
%   yield the same numerical result. The algorithms differ in run-time
%   efficiency, as well as in the maximum pattern order ORD supported.
%
%   Possible values for ALGORITHM are:
%
%     ---------------------------------------------------------------------
%     ALGORITHM    MAX ORD    NEEDS MEX    REMARKS
%     ---------------------------------------------------------------------
%     lookup            10           no
%     overlap           20           no
%     plain             20           no    Very slow.
%     vectorised        20           no    Best choice if MEX unavailable.
%     ---------------------------------------------------------------------
%     lookup_c          10          yes    Fast for small order ORD.
%     overlap_c         20          yes    Good general-purpose algorithm.
%     overlap_mp_c     255          yes    Supports orders > 20.
%     plain_c           20          yes
%     ---------------------------------------------------------------------
%
%   The dimensionality of the returned array N is equivalent to X (except
%   for the 'overlap_mp_c' algorithm, see below). The ordinal patterns
%   in N are encoded by integers ranging from 0 to factorial(ORD)-1.
%   The data type of N depends on the pattern order ORD:
%
%     For  1 < ORD <  19, 'double' values are used.
%     For 18 < ORD <  21, 'uint64' values are used.
%     For 20 < ORD < 255, 'uint64' arrays are used.
%
%   The 'overlap_mp_c' algorithm supports pattern orders ORD up to 255.
%   When using this algorithm, the first dimension of the array N indexes
%   the 64-bit wide elements of single ordinal pattern. For example,
%
%     >>  'ordpat(zeros(1000, 8), 21, 1, 'overlap_mp_c')'
%
%   returns a uint64-valued array of size [2, 980, 8], because a sequence
%   of length 1000 has 980 ordinal patterns of order ORD == 21, whereas
%   storing each such pattern requires two 64-bit values.
%
%   See the publication: Berger S, Kravtsiv A, Schneider G, Jordan D.
%                        Teaching Ordinal Patterns to a Computer.
%                        Entropy. 2019; 21(10):1023.
%                        https://doi.org/10.3390/e21101023

%   Copyright (c) 2019, Sebastian Berger.
%
%   Klinikum rechts der Isar der
%   Technischen Universitaet Muenchen
%   Munich, Germany
%
%   All rights reserved.
%
%   Redistribution and use in source and binary forms, with or without
%   modification, are permitted provided that the following conditions are
%   met:
%       * Redistributions of source code must retain the above copyright
%         notice, this list of conditions and the following disclaimer.
%       * Redistributions in binary form must reproduce the above copyright
%         notice, this list of conditions and the following disclaimer in
%         the documentation and/or other materials provided with the
%         distribution.
%       * Neither the names of the copyright holders nor the names of its
%         contributors may be used to endorse or promote products derived
%         from this software without specific prior written permission.
%
%   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
%   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
%   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
%   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR OR
%   THE KLINIKUM RECHTS DER ISAR BE LIABLE FOR ANY DIRECT, INDIRECT,
%   INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
%   BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
%   LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
%   AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
%   OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
%   THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
%   DAMAGE.


% for storing whether MEX functions are available
persistent has_mex;

% for storing the most recently used lookup table
persistent lookup;

if nargin == 1
    if strcmp(x, '--make')
        make();
        has_mex = is_mex_available();
        return;
    elseif strcmp(x, '--clean')
        make('clean');
        has_mex = [];
        return;
    end
end % fall through if neither '--make' nor '--clean'


narginchk(2, 4);

if isempty(has_mex)
    has_mex = is_mex_available();
end

% default algorithms
if nargin < 4
    if ord <= 20
        if has_mex
            algorithm = 'overlap_c';
        else
            algorithm = 'vectorised';
        end
    else
        if has_mex
            algorithm = 'overlap_mp_c';
        else
            error(['ORD > 20 requires C library functions. ', ...
                   'Please call ''ordpat --make''.']);
        end
    end
end

% default time lag
if nargin < 3 || isempty(lag)
    lag = 1;
end

% throws error if arguments invalid
check_input_arguments(x, ord, lag, algorithm, has_mex);

% convert to 'double' because the C library does not support other formats
x = double(x);

% select the function handle
enc_function = str2func(['encode_', algorithm]);

% prepare lookup table if required
if ismember(algorithm, {'lookup', 'lookup_c'})
    lookup = update_lookup_table(lookup, ord);
    enc_function = @(x, ord, lag) enc_function(x, ord, lag, lookup.tab);
end

% Only the 'encode_vectorised' function supports multi-dimensional input.
% In all other cases: (1) interpret 'x' as a vector,
%                     (2) encode it,
%                     (3) truncate the edges, and
%                     (4) reshape the result.
switch algorithm
    case 'vectorised'
        n = enc_function(x, ord, lag);

    otherwise
        multi_precision = strcmp(algorithm, 'overlap_mp_c');

        if multi_precision
            width = pattern_width(ord);
        else
            width = 1;
        end

        % initialise (slightly oversized) results matrix
        dims    = size(x);
        dims(1) = dims(1) * width;

        if multi_precision || ord > 18
            n = uint64(zeros(dims));
        else
            n = zeros(dims);
        end

        % number of patterns if 'x' would be a vector
        num_pat = width * (numel(x) - (ord - 1) * lag);

        % extract and encode ordinal patterns of vector 'x(:)'
        n(1:num_pat) = enc_function(x(:), ord, lag);

        % truncate
        dims    = size(x);
        dims(1) = dims(1) - (ord - 1) * lag;
        num_pat = dims(1) * width;

        n = n(1:num_pat, :);

        % reshape
        if multi_precision
            shape = [width, dims];
        else
            shape = dims;
        end

        n = reshape(n, shape);
end


%%% Make sure the input arguments are valid, otherwise throw an error.
function check_input_arguments(x, ord, lag, algorithm, has_mex)

% the algorithms, and their supported pattern orders
algorithms = struct('lookup',       {{ 10, false}}, ...
                    'lookup_c',     {{ 10,  true}}, ...
                    'overlap',      {{ 20, false}}, ...
                    'overlap_c',    {{ 20,  true}}, ...
                    'overlap_mp_c', {{255,  true}}, ...
                    'plain',        {{ 20, false}}, ...
                    'plain_c',      {{ 20,  true}}, ...
                    'vectorised',   {{ 20, false}});

if ~isreal(ord) || ~isscalar(ord) || ord < 2 || floor(ord) ~= ord
    error('ORD must be an integer greater than 1.');
end

if ~isreal(lag) || ~isscalar(lag) || lag < 1 || floor(lag) ~= lag
    error('LAG must be a positive integer.');
end

if ~isnumeric(x) || ~all(isfinite(x(:))) || any(imag(x(:)))
    error('X must contain real, finite, numerical data.');
end

if (size(x, 1) - (ord - 1) * lag) < 1
    error('Input data X too short for selected ORD and LAG.');
end

if ~ischar(algorithm) || ~isrow(algorithm)
    error('ALGORITHM must be a character string.');
end

% check if algorithm is known and get its maximum pattern order
try
    info = algorithms.(algorithm);
catch
    error('Unknown algorithm: ''%s''.', algorithm);
end

if ord > info{1}
    error(['The ''%s'' algorithm does not support orders ' ...
           'greater than %d.'], algorithm, info{1});
end

if ~has_mex && info{2}
    error(['The ''%s'' algorithm requires C library functions. ', ...
           'Please call ''ordpat --make''.'], algorithm);
end

% make sure we are not loosing precision
if any(abs(x(:)) > flintmax())
    error('X cannot be converted to ''double'' without loss of precision.')
end


%%% Check if MEX functions are available.
function tf = is_mex_available()

tf = true;
try
    pattern_width(2);
catch err
    tf = false;
end


%%% Recreate the lookup table if necessary.
function lookup = update_lookup_table(lookup, ord)

if isempty(lookup) || lookup.ord ~= ord
    try
        lookup.tab = create_lookup_table_c(ord);
    catch
        lookup.tab = create_lookup_table(ord);
    end

    lookup.ord = ord;
end
