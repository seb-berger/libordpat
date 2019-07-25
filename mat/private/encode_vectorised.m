function y = encode_vectorised(x, ord, lag)
%ENCODE_VECTORISED Implementation of the 'plain' encoding algorithm.
%   N = ENCODE_VECTORISED(X, ORD, LAG) extracts the consecutive ordinal
%   patterns of order ORD from the time series contained in X, using the
%   time lag LAG and a vectorised version of the 'plain' encoding
%   algorithm.
%
%   Orders 1 < ORD < 21, and any time lag LAG > 0 are supported.
%
%   ENCODE_VECTORISED operates along the first dimension of X, that is,
%   each column of X is treated as an individual time series.
%
%   The dimensionality of the returned array N is equivalent to X.
%   Its data format is 'double' if ORD < 19, and it is 'uint64' for
%   18 < ORD < 21. The ordinal patterns in N are encoded by integers
%   ranging from 0 to factorial(ORD) - 1.
%
%   WARNING: THIS IS AN INTERNAL FUNCTION WITH MINIMAL ERROR CHECKING!

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

persistent ctx;

dims = size(x);
dims(1) = dims(1) - (ord - 1) * lag;

if isempty(ctx) || ctx.ord ~= ord || ctx.lag ~= lag

    % create embedding code
    cmd = '';
    for i = 1:ord
        i1  = (i   - 1) * lag + 1;
        i2  = (ord - i) * lag;
        cmd = sprintf('%sx%d = x(%d : end - %d, :);\n', cmd, i, i1, i2);
    end

    % create Lehmer code
    cmd = sprintf('%s\ny = ', cmd);

    if ord <= 18
        format = '%s%d * ((x%d > x%d)';
    else
        format = '%s%d * uint64((x%d > x%d)';
    end

    for i = 1:(ord - 1)
        fct = factorial(ord - i);
        cmd = sprintf(format, cmd, fct, i, i + 1);

        for j = (i + 2):ord
            cmd = sprintf('%s + (x%d > x%d)', cmd, i, j);
        end

        if i == ord - 1
            cmd = sprintf('%s);\n', cmd);
        else
            cmd = sprintf('%s) + ...\n%s', cmd, repmat(' ', 1, 4));
        end
    end

    ctx.ord = ord;
    ctx.lag = lag;
    ctx.cmd = cmd;
end

% execute the code
eval(ctx.cmd);
y = reshape(y, dims);
