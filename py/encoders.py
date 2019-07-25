#   Copyright (c) 2019, Sebastian Berger.
#
#   Klinikum rechts der Isar der
#   Technischen Universitaet Muenchen
#   Munich, Germany
#
#   All rights reserved.
#
#   Redistribution and use in source and binary forms, with or without
#   modification, are permitted provided that the following conditions are
#   met:
#       * Redistributions of source code must retain the above copyright
#         notice, this list of conditions and the following disclaimer.
#       * Redistributions in binary form must reproduce the above copyright
#         notice, this list of conditions and the following disclaimer in
#         the documentation and/or other materials provided with the
#         distribution.
#       * Neither the names of the copyright holders nor the names of its
#         contributors may be used to endorse or promote products derived
#         from this software without specific prior written permission.
#
#   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
#   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR OR
#   THE KLINIKUM RECHTS DER ISAR BE LIABLE FOR ANY DIRECT, INDIRECT,
#   INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
#   BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
#   LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
#   AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
#   OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
#   THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
#   DAMAGE.
"""
libordpat - Native implementations of ordinal pattern encoding algorithms.

This module provides Python implementations of the

  * 'plain'
  * 'overlap'
  * 'lookup'
  * 'vectorised'

algorithms, which all allow for encoding ordinal patterns of order m
into non-negative integers between 0 and factorial(ord) - 1.

The 'vectorised' algorithm is probably the only one that makes practical
sense for a Python implementation.

"""
import numpy as _np
from math import factorial as _factorial

import utils as _utils


def _encode_pattern(x):
    ord = len(x)      # pattern order
    n_rel = ord - 1   # number of order relations

    y = 0

    for i in range(0, n_rel):
        for j in range(i + 1, ord):
            y += x[i] > x[j]

        y *= n_rel - i

    return y


def create_lookup_table(ord):
    """
    Return pattern table to be used with the 'lookup' algorithm.

    The pattern order must be an integer greater than 1. Because the
    table size increases super-exponentially with the pattern order,
    only orders `ord` < 11 are accepted. (A lookup table for `ord` == 11
    would occupy more than 3 GiB of memory).

    Parameters
    ----------
    ord : int
        Pattern order between 2 and 10.

    Returns
    -------
    table : ndarray of uint64
        A lookup table containing factorial(ord) * ord pattern codes.
        Each code is a value between 0 and factorial(ord) - 1.

    """
    if not _utils.is_scalar_int(ord):
        raise TypeError("order must be a scalar integer")

    if ord < 2 or ord > 10:
        raise ValueError("order must be between 2 and 10")

    pat = 2 * _utils.pattern_table(ord - 1)
    pat = pat.repeat(ord, 0)

    val = _np.arange(2.0*ord - 1, 0, -2)
    val = val.reshape((ord, 1))
    val = _np.tile(val, (_factorial(ord -1), 1))

    pat = _np.concatenate((pat, val), axis=1)

    tab = encode_vectorised(pat, ord, 1, trustme=True)
    tab = tab.reshape((-1, ord))
    tab = _np.tile(tab, (ord, 1))

    return tab


def encode_plain(x, ord, lag, trustme=False):
    """
    Extract and encode ordinal patterns using the 'plain' algorithm.

    Parameters
    ----------
    x : array_like
        Time series data to be encoded. Can be multi-dimensional, and must be
        convertible to an ndarray of data type 'double'. The input array is
        flattened prior to encoding.

    ord : int
        Order of the encoding, an integer between 2 and 20.

    lag : int
        Time lag used for creating embedding vectors. Must be a positive
        integer.

    trustme : bool, optional
        If set to True, the function arguments are not validated before
        processing. This will speed up the calculation.

    Returns
    -------
    pat : ndarray of uint64
        One-dimensional array, with each value representing an ordinal
        pattern. It holds that ``pat.size == x.size - (ord - 1) * lag``.

    """
    if not trustme:
        x = _utils.prepare_input_data(x, ord, lag)
        if ord > 20:
            raise ValueError("plain algorithm does not support ord > 20")

    n_rel = ord - 1                # number of order relations
    n_pat = x.size - n_rel * lag   # sequence length

    pat = _np.zeros(n_pat, dtype=_np.uint64)

    for k in range (0, n_pat):
        for i in range(0, n_rel):
            for j in range(i + 1, ord):
                pat[k] += x[k + i*lag] > x[k + j*lag]

            pat[k] *= n_rel - i

    return pat


def encode_overlap(x, ord, lag, trustme=False):
    """
    Extract and encode ordinal patterns using the 'overlap' algorithm.

    Parameters
    ----------
    x : array_like
        Time series data to be encoded. Can be multi-dimensional, and must be
        convertible to an ndarray of data type 'double'. The input array is
        flattened prior to encoding.

    ord : int
        Order of the encoding, an integer between 2 and 20.

    lag : int
        Time lag used for creating embedding vectors. Must be a positive
        integer.

    trustme : bool, optional
        If set to True, the function arguments are not validated before
        processing. This will speed up the calculation.

    Returns
    -------
    pat : ndarray of uint64
        One-dimensional array, with each value representing an ordinal
        pattern. It holds that ``pat.size == x.size - (ord - 1) * lag``.

    """
    if not trustme:
        x = _utils.prepare_input_data(x, ord, lag)
        if ord > 20:
            raise ValueError("overlap algorithm does not support ord > 20")

    n_rel = ord - 1                # number of order relations
    n_pat = x.size - n_rel * lag   # sequence length

    pat = _np.zeros(n_pat, dtype=_np.uint64)

    ranks = _np.zeros([lag, n_rel], dtype=_np.uint8)

    radix = [_factorial(r) for r in range(n_rel, 0, -1)]
    radix = _np.array(radix, dtype=_np.uint64)

    for k in range(0, lag):
        for i in range(0, n_rel - 1):
            for j in range(i + 1, n_rel):
                ranks[k, i + 1] += x[k + i*lag] > x[k + j*lag]

    i = 0

    for k in range(0, n_pat):
        ranks[i] = _np.roll(ranks[i], -1)
        ranks[i, -1] = 0

        updates = x[k:k + n_rel*lag:lag] > x[k + n_rel*lag]
        ranks[i] += updates

        pat[k] = _np.matmul(radix, ranks[i])
        i = (i + 1) % lag

    return pat


def encode_lookup(x, ord, lag, tab, trustme=False):
    """
    Extract and encode ordinal patterns using the 'lookup' algorithm.

    Parameters
    ----------
    x : array_like
        Time series data to be encoded. Can be multi-dimensional, and must be
        convertible to an ndarray of data type 'double'. The input array is
        flattened prior to encoding.

    ord : int
        Order of the encoding, an integer between 2 and 10.

    lag : int
        Time lag used for creating embedding vectors. Must be a positive
        integer.

    tab : ndarray of uint64
        A lookup table containing factorial(ord) * ord pattern codes, as
        returned by the `create_lookup_table` function. Note that custom
        mappings are NOT supported, and will lead to undefined behaviour.

    trustme : bool, optional
        If set to True, the function arguments are not validated before
        processing. This will speed up the calculation.

    Returns
    -------
    pat : ndarray of uint64
        One-dimensional array, with each value representing an ordinal
        pattern. It holds that ``pat.size == x.size - (ord - 1) * lag``.

    """
    if not trustme:
        x = _utils.prepare_input_data(x, ord, lag)
        if ord > 10:
            raise ValueError("lookup algorithm does not support ord > 10")

        _utils.check_lookup_table(tab, ord)

    n_rel = ord - 1
    n_pat = x.size - n_rel * lag

    pat = _np.zeros(n_pat, dtype=_np.uint64)

    for i in range(0, lag):
        pat[i] = _encode_pattern(x[i:i + ord*lag:lag])

    next_val = x[n_rel*lag:]
    ranks = _np.zeros(next_val.shape, _np.uint8)

    for i in range(0, n_rel):
        ranks += x[i*lag : (i - n_rel)*lag] > next_val

    for i in range(lag, n_pat):
        pat[i] = tab[pat[i - lag], ranks[i]]

    return pat


_scope = {'order': None, 'lag': None, 'cmd': None}


def encode_vectorised(x, ord, lag, axis=-1, trustme=False):
    """
    Extract and encode ordinal patterns using the 'vectorised' algorithm.

    Parameters
    ----------
    x : array_like
        Time series data to be encoded. Can be multi-dimensional, and must be
        convertible to an ndarray of data type 'double'.

    ord : int
        Order of the encoding, an integer between 2 and 20.

    lag : int
        Time lag used for creating embedding vectors. Must be a positive
        integer.

    axis : int, optional
        The axis of the input data along which the extraction shall be
        performed. By default, the last axis is used.

    trustme : bool, optional
        If set to True, the function arguments are not validated before
        processing. This will speed up the calculation.

    Returns
    -------
    pat : ndarray of uint64
        Each value represents an ordinal pattern.  The shape of the `pat`
        array is equal to the input array `x`, except for the axis on
        which the encoding is performed. For this axis, it holds that

          ``pat.shape[axis] == x.shape[axis] - (ord - 1) * lag``.

    """
    if not trustme:
        x = _utils.prepare_input_data(x, ord, lag, axis)
        if ord > 20:
            raise ValueError("vectorised algorithm does not support ord > 20")

    _scope['x'] = _np.moveaxis(x, axis, 0)

    if (_scope['order'] != ord) or (_scope['lag'] != lag):
        cmd = ""
        for idx in range(1, ord + 1):
            fwd = (idx - 1)*lag;
            rev = (ord - idx)*lag;
            cmd += "x{0} = x[{1}:{2}]\n".format(idx, fwd, -rev or None)

        cmd += "\ny = "

        word_size = _utils.pattern_word_size(ord)
        for i in range(1, ord):
            cmd += "{0} * (".format(_factorial(ord - i))
            cmd += "(x{0} > x{1})".format(i, i + 1)
            cmd += ".astype(_np.uint8)"

            for j in range(i + 2, ord + 1):
                cmd += " + (x{0} > x{1})".format(i, j)

            cmd += ").astype(_np.uint{0})".format(word_size)

            if i == ord - 1:
                cmd += "\n"
            else:
                cmd += " + \\\n" + " "*4

        _scope['cmd'] = compile(cmd, "<string>", "exec")
        _scope['order'] = ord
        _scope['lag'] = lag

    exec(_scope['cmd'], None, _scope)

    pat = _scope['y']
    pat = _np.moveaxis(pat, 0, axis)

    return pat.astype(_np.uint64)
