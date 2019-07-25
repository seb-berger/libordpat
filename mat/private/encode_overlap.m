function n = encode_overlap(x, ord, lag)
%ENCODE_OVERLAP Implementation of the 'overlap' encoding algorithm.
%   N = ENCODE_OVERLAP(X, ORD, LAG) extracts the consecutive ordinal
%   patterns of order ORD from the time series contained in X, using the
%   time lag LAG and the 'overlap' encoding algorithm.
%
%   Orders 1 < ORD < 21, and any time lag LAG > 0 are supported.
%
%   ENCODE_OVERLAP operates along a flattend version of X, that is,
%   the whole array X is treated as one time series.
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

N = numel(x) - (ord - 1) * lag;
r = zeros(lag, ord);

if ord <= 18
    n = zeros(N, 1);
else
    n = uint64(zeros(N, 1));
end

% extract inversion counts for the first LAG patterns
for k = 1:lag
    for i = 1:(ord - 2)
        for j = (i + 1):(ord - 1)
            r(k, i + 1) = r(k, i + 1) + (x(k + (i - 1) * lag) > ...
                                         x(k + (j - 1) * lag));
        end
    end
end

% encode X
i = 1;
for k = 1:N
    for j = 1:(ord - 1)
        % update inversion counts
        r(i, j) = r(i, j + 1) + (x(k + (j - 1)   * lag) > ...
                                 x(k + (ord - 1) * lag));

        % encode the pattern
        n(k) = (ord - j) * (n(k) + r(i, j));
    end

    i = mod(i, lag) + 1;
end
