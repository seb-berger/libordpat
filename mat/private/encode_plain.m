function n = encode_plain(x, ord, lag)
%ENCODE_PLAIN Implementation of the 'plain' encoding algorithm.
%   N = ENCODE_PLAIN(X, ORD, LAG) extracts the consecutive ordinal
%   patterns of order ORD from the time series contained in X, using the
%   time lag LAG and the 'plain' encoding algorithm.
%
%   Orders 1 < ORD < 21, and any time lag LAG > 0 are supported.
%
%   ENCODE_PLAIN operates along a flattend version of X, that is, the
%   whole array X is treated as one time series.
%
%   The returned value N is a column vector. Its data format is 'double'
%   if ORD < 19, and it is 'uint64' for 18 < ORD < 21. The ordinal patterns
%   in N are encoded by integers ranging from 0 to factorial(ORD) - 1.
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

if ord <= 18
    n = enc_double(x(:), ord, lag);
else
    n = enc_uint64(x(:), ord, lag);
end


%%% Encode patterns using a 'double' output array.
function n = enc_double(x, ord, lag)

N = numel(x) - (ord - 1) * lag;
n = zeros(N, 1);

for k = 1:N
    for i = 1:(ord - 1)
        for j = (i + 1):ord
            n(k) = n(k) + (x(k + (i - 1) * lag) > x(k + (j - 1) * lag));
        end

        n(k) = n(k) * (ord - i);
    end
end


%%% Encode patterns using a 'uint64' output array.
function n = enc_uint64(x, ord, lag)

N = numel(x) - (ord - 1) * lag;
n = uint64(zeros(N, 1));

for k = 1:N
    for i = 1:(ord - 1)
        for j = (i + 1):ord
            n(k) = n(k) + uint64(x(k + (i - 1) * lag) > ...
                                 x(k + (j - 1) * lag));
        end

        n(k) = n(k) * (ord - i);
    end
end
