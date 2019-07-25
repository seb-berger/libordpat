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
libordpat - Utilities for the libordpat library.

This module is a small collection of nifty support functions for the
ordpat module. Some may also be useful in their own right, though.

"""
from __future__ import division as _division

import numpy as _np
from math import factorial as _factorial, ceil as _ceil, log as _log
from itertools import permutations as _permutations


def is_scalar_int(x):
    """
    Return True if `x` represents a scalar integer value.

    Parameters
    ----------
    x : array_like
        Argument to be tested.

    Returns
    -------
    tf : bool
        True if `x` is a scalar integer, False otherwise.

    """
    return _np.ndim(x) == 0 and int(x) == x


def check_lookup_table(tab, ord):
    """
    Raise an error if `tab` is not a valid lookup table.

    Parameters
    ----------
    tab : possibly a duck
        Python object to be tested.

    ord : int
        Pattern order between 2 and 10.

    """
    if not is_scalar_int(ord):
        raise TypeError("ord must be a scalar integer")

    if ord < 2 or ord > 10:
        raise ValueError("ord must be between 2 and 10")

    if not isinstance(tab, _np.ndarray):
        raise TypeError("tab must be an ndarray")

    if not tab.dtype == _np.uint64:
        raise TypeError("tab must have dtype uint64")

    if not tab.flags["CA"]:
        raise TypeError("tab must be contiguous and word-aligned")

    if not tab.shape == (_factorial(ord), ord):
        raise ValueError("tab must have shape ord! x ord")

    # Invalid codes may cause segmentation fault.
    if (tab >= tab.shape[0]).any():
        raise ValueError("tab contains invalid pattern codes")


def prepare_input_data(x, ord, lag, axis=None):
    """
    Prepare a set of input arguments for encoding.

    Checks input arguments for validity/consistency. If successful,
    converts `x` into an ndarray of data type double.

    Parameters
    ----------
    x : array_like
        Input data.

    ord : int
        Pattern order between 2 and 255.

    lag : int
        Time lag, greater than 0.

    axis : int, optional
        Axis along which the encoding will be performed.

    Returns
    -------
    x : ndarray of double
        Array containing the data in `x`, stored in row-major (C style)
        order.  If `axis` is None, then a flattened 1D array is returned,
        otherwise the shape of the input data is preserved.

    """
    if (axis is not None) and (not is_scalar_int(axis)):
        raise TypeError("axis must be None or a scalar integer")

    if not is_scalar_int(ord):
        raise TypeError("order must be a scalar integer")

    if ord < 2 or ord > 255:
        raise ValueError("order must be between 2 and 255")

    if not is_scalar_int(lag):
        raise TypeError("time lag must be a scalar integer")

    if lag < 1:
        raise ValueError("time lag must be positive")

    x = _np.require(x, _np.double, "CA")

    if axis is None:
        x = x.reshape(-1)
        axis = 0

    if x.shape[axis] <= (ord - 1) * lag:
        raise ValueError("x is too short for given order and time lag")

    return x


def pattern_word_size(ord):
    """
    Return word size of an ordinal pattern of order ord.

    Depending on the order, the numerical representation of an ordinal pattern
    requires a certain number of bits. This function rounds the bit width up
    to the next-greatest machine word.

    Parameters
    ----------
    ord : int
        Pattern order between 2 and 20.

    Returns
    -------
    width : int
        Either 8, 16, 32 or 64.

    """
    if not is_scalar_int(ord):
        raise TypeError("order must be a scalar integer")

    if ord < 2 or ord > 20:
        raise ValueError("order must be between 2 and 20")

    n_pats = _factorial(ord)
    n_bits = _log(n_pats, 2)
    exp    = _ceil(_log(n_bits, 2))

    return int(max(8, 2**exp))


def pattern_uint64_width(ord):
    """
    Return number of 64-bit words occupied by a pattern of order `ord`.

    Depending on the order, the numerical representation of an ordinal pattern
    requires a certain number of bits. This function rounds the bit width up
    to the next-greatest multiple of a 64-bit machine word.

    Parameters
    ----------
    ord : int
        Pattern order between 2 and 255.

    Returns
    -------
    width : int
        An integer between 1 and 27.

    """
    if not is_scalar_int(ord):
        raise TypeError("order must be a scalar integer")

    if ord < 2 or ord > 255:
        raise ValueError("order must be between 2 and 255")

    n_pats  = _factorial(ord)
    n_bits  = _log(n_pats, 2)
    n_words = _ceil(n_bits / 64)

    return int(n_words)


def pattern_table(ord):
    """
    Return all ordinal patterns of order `ord` in their rank representation.

    Parameters
    ----------
    ord : int
        Pattern order between 2 and 10.

    Returns
    -------
    tab : ndarray of double
        2D array, with each row representing an ordinal pattern in rank its
        representation. Patterns are sorted in lexicographic order.

    """
    if not is_scalar_int(ord):
        raise TypeError("order must be a scalar integer")

    if ord < 1 or ord > 10:
        raise ValueError("order must be between 2 and 10")

    tmp = range(1, ord + 1)
    return _np.asarray([i for i in _permutations(tmp)], dtype=_np.double)
