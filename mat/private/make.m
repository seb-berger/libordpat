function make(target)
%MAKE Build C library functions.
%   MAKE(TARGET) prepares the C library functions required by ordpat.
%   It uses the MEX build architecture. TARGET can either be:
%
%     - 'all'   to create the MEX files needed, or
%     - 'clean' to delete previously created MEX files.
%
%   MAKE without an argument is equivalent to MAKE('all').
%
%   If the MAKE function does not succeed, ensure that a compatible
%   C compiler is installed, and that MEX is configured correctly.

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

if nargin == 0
    target = 'all';
end

args.optimisation = '-O3 -std=c99 -Wall -Wextra -pedantic';
args.lib_path     = abspath(fullfile('..', '..', 'lib'));
args.src_path     = abspath(fullfile('..',  'mex'));
args.out_path     = abspath(fullfile('.'));
args.include      = ['-I', args.lib_path];

switch target
    case 'all'
        make_all(args);
    case 'clean'
        make_clean(args);
    otherwise
        error('Unknown target: %s.', target);
end


%%% Compile all MEX files.
function make_all(args)

opt = args.optimisation;
inc = args.include;
src = args.src_path;
lib = args.lib_path;
out = args.out_path;

if exist('OCTAVE_VERSION', 'builtin') ~= 0
    mex(opt,                             ... % optimisation
        '-c',                            ... % compile only
        '-DBUILD_LIBORDPAT_STATIC',      ... % no export for ordpat.o
        '-o', fullfile(out, 'ordpat.o'), ... % output file
        fullfile(lib, 'ordpat.c'));          % source file

    build = @(dst, file)                      ...
            mex(opt,                          ... % optimisation
                '-DBUILD_LIBORDPAT_STATIC',   ... % no export for ordpat.o
                '-s',                         ... % strip binary
                inc,                          ... % includes
                '-o', fullfile(out, dst),     ... % target dir
                fullfile(out, 'ordpat.o'),    ... % ordpat.o obj
                fullfile(src, file));             % source file

    build('create_lookup_table_c', 'create_lookup_table.c');
    build('encode_plain_c',        'encode_plain.c');
    build('encode_overlap_c',      'encode_overlap.c');
    build('encode_overlap_mp_c',   'encode_overlap_mp.c');
    build('encode_lookup_c',       'encode_lookup.c');
    build('pattern_width',         'pattern_width.c');
    build('xorshift',              'xorshift.c');

    delete(fullfile(out, 'ordpat.o'));

else % MATLAB
    mex(['COPTIMFLAGS="', opt, ' -DNDEBUG"'], ... % optimisation
        '-DBUILD_LIBORDPAT_STATIC',           ... % no export for ordpat.o
        '-silent',                            ... % suppress prompt
        '-largeArrayDims',                    ... % 64-bit mwSize
        '-c',                                 ... % compile only
        '-outdir', out,                       ... % target dir
        fullfile(lib, 'ordpat.c'));               % source file

    if ispc()
        obj_file = 'ordpat.obj';
    else
        obj_file = 'ordpat.o';
    end

    build = @(dst, file)                              ...
            mex(['COPTIMFLAGS="', opt, ' -DNDEBUG"'], ... % optimisation
                '-DBUILD_LIBORDPAT_STATIC',           ... % no export for ordpat.o
                '-silent',                            ... % suppress prompt
                '-largeArrayDims',                    ... % 64-bit mwSize
                inc,                                  ... % includes
                '-outdir', out,                       ... % target dir
                '-output', dst,                       ... % output file
                fullfile(out, obj_file),              ... % ordpat.o
                fullfile(src, file));                     % source file

    build('create_lookup_table_c', 'create_lookup_table.c');
    build('encode_plain_c',        'encode_plain.c');
    build('encode_overlap_c',      'encode_overlap.c');
    build('encode_overlap_mp_c',   'encode_overlap_mp.c');
    build('encode_lookup_c',       'encode_lookup.c');
    build('pattern_width',         'pattern_width.c');
    build('xorshift',              'xorshift.c');

    delete(fullfile(out, obj_file));
end


%%% Delete all MEX files.
function make_clean(args)

clear('-f'); % unload

for suffix = {mexext(), 'o', 'obj'}
    files = dir(fullfile(args.out_path, ['*.', suffix{:}]));
    files = fullfile(args.out_path, {files.name});
    cellfun(@delete, files);
end


%%% Return absolute path.
function absolute = abspath(relative)

origin = fileparts(mfilename('fullpath'));

absolute = what(fullfile(origin, relative));
if isempty(absolute)
    error('Invalid path: %s', relative);
end

absolute = absolute.path;
