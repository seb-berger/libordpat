function tab = create_lookup_table(ord)
%CREATE_LOOKUP_TABLE Table for the 'lookup' and 'lookup_c' algorithms.
%   TAB = CREATE_LOOKUP_TABLE(ORD) creates a lookup table TAB for
%   encoding ordinal patterns of order ORD using either the 'lookup'
%   or 'lookup_c' pattern encoding algorithms.
%
%   The pattern order ORD must be an integer greater than 1. Because the
%   table size increases super-exponentially with the pattern order,
%   only orders ORD < 11 are accepted. (A lookup table for ORD == 11
%   would occupy more than 3 GiB of memory).

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

if ~isscalar(ord) || floor(ord) ~= ord || ~isreal(ord) || ord < 2
    error('ORD must be a scalar integer greater than 1.');
end

if ord > 10
    error('ORD is limited to a maximum of 10.');
end

% This is just a bunch of combinatorial tricks and rearrangements, easily
% comprehended by stepping through the code with a debugger.
pat = 2 * flipud(perms(1:ord - 1));
pat = kron(pat, ones(ord, 1));

val = (2 * ord - 1):-2:1;
val = repmat(val.', factorial(ord - 1), 1);

tab = encode_vectorised([pat, val].', ord, 1);
tab = reshape(tab, ord, []);
tab = repmat(tab, 1, ord);
tab = uint64(tab);
