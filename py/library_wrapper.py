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
libordpat - Python wrapper around the libordpat library.

This module wraps the ordinal pattern encoding functions of the
'libordpat' C library for use in Python.

Upon import of the module, the library file

  * libordpat.dll    (Windows)
  * libordpat.dylib  (Mac OS X)
  * libordpat.so     (Linux and other Unices)

must be installed somewhere in the system's library search path.

"""
from __future__ import division as _division

import os as _os
import sys as _sys
import numpy as _np
import ctypes as _ctypes
from math import factorial as _factorial

import utils as _utils

_mod_path = _os.path.dirname(__file__)
_filename = _os.path.join(_mod_path, "..", "lib", "libordpat")

if _os.name == "nt":
    _filename += ".dll"
elif _os.name == "posix":
    if _sys.platform.startswith("darwin"):
        _filename += ".dylib"
    else:
        _filename += ".so.0.1.0"
else:
    raise OSError("operating system not supported")

# two useful ndpointer instances
_ptr_uint64 = _np.ctypeslib.ndpointer(dtype=_np.uint64, flags=("C", "A"))
_ptr_double = _np.ctypeslib.ndpointer(dtype=_np.double, flags=("C", "A"))

# try to load the shared library
_lib = _ctypes.cdll[_filename]

_lib.ordpat_create_lookup_table.argtypes = [
    _ptr_uint64,      # pointer to output array
    _ctypes.c_size_t, # size of array
    _ctypes.c_uint    # pattern order
    ]

_lib.ordpat_xorshift_rand_double.restype  = None
_lib.ordpat_xorshift_rand_double.argtypes = [
        _ptr_double,      # pointer to output array
        _ctypes.c_size_t, # size of array
        _ctypes.c_uint32  # seed
    ]

# Most encoding functions have the same arguments:
_encode_argtypes = [
    _ptr_double,      # input array
    _ctypes.c_size_t, # size of input array
    _ptr_uint64,      # output array
    _ctypes.c_size_t, # size of output array
    _ctypes.c_uint,   # pattern order
    _ctypes.c_uint    # time lag
    ]

_lib.ordpat_encode_plain.argtypes      = _encode_argtypes
_lib.ordpat_encode_overlap.argtypes    = _encode_argtypes
_lib.ordpat_encode_overlap_mp.argtypes = _encode_argtypes

# In addition, 'ordpat_encode_lookup' requires a pointer to a lookup table.
_lib.ordpat_encode_lookup.argtypes = _encode_argtypes + [_ptr_uint64]


def _raise_api_error(ret):
    msg = ("ORDPAT_SUCCESS",
           "ORDPAT_ERROR_ORDER_INVALID",
           "ORDPAT_ERROR_LAG_INVALID",
           "ORDPAT_ERROR_INPUT_TOO_SHORT",
           "ORDPAT_ERROR_OUTPUT_TOO_SHORT")

    raise RuntimeError("library call returned: " + msg[ret])


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

    table = _np.empty((_factorial(ord), ord), dtype=_np.uint64)
    ret = _lib.ordpat_create_lookup_table(table, table.size, ord)

    if not ret == 0:
        _raise_api_error(ret)

    return table


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

    n_pat = x.size - (ord - 1) * lag

    pat = _np.empty(n_pat , dtype=_np.uint64)

    ret = _lib.ordpat_encode_plain(x, x.size, pat, pat.size, ord, lag)

    if not ret == 0:
        _raise_api_error(ret)

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

    n_pat = x.size - (ord - 1) * lag

    pat = _np.empty(n_pat , dtype=_np.uint64)

    ret = _lib.ordpat_encode_overlap(x, x.size, pat, pat.size, ord, lag)

    if not ret == 0:
        _raise_api_error(ret)

    return pat


def encode_overlap_mp(x, ord, lag, trustme=False):
    """
    Extract and encode ordinal patterns using the 'overlap_mp' algorithm.

    Parameters
    ----------
    x : array_like
        Time series data to be encoded. Can be multi-dimensional, and must be
        convertible to an ndarray of data type 'double'. The input array is
        flattened prior to encoding.

    ord : int
        Order of the encoding, an integer between 2 and 255.

    lag : int
        Time lag used for creating embedding vectors. Must be a positive
        integer.

    trustme : bool, optional
        If set to True, the function arguments are not validated before
        processing. This will speed up the calculation.

    Returns
    -------
    pat : ndarray of uint64
        Two-dimensional array, with each vector along the second axis
        representing an ordinal pattern. It holds that

          ``pat.shape[0] == x.size - (ord - 1) * lag``,

        whereas ``pat.shape[1]`` depends on the pattern order.

    """
    if not trustme:
        x = _utils.prepare_input_data(x, ord, lag)

        if ord > 255:
            raise ValueError("overlap_mp algorithm does not support ord > 255")

    n_pat = x.size - (ord - 1) * lag

    width = _utils.pattern_uint64_width(ord)
    n_pat *= width

    pat = _np.empty(n_pat, dtype=_np.uint64)

    ret = _lib.ordpat_encode_overlap_mp(x, x.size, pat, pat.size, ord, lag)

    if not ret == 0:
        _raise_api_error(ret)

    pat = pat.reshape((-1, width))

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

    n_pat = x.size - (ord - 1) * lag

    pat = _np.empty(n_pat , dtype=_np.uint64)

    ret = _lib.ordpat_encode_lookup(x, x.size, pat, pat.size, ord, lag, tab)

    if not ret == 0:
        _raise_api_error(ret)

    return pat


def xorshift_rand(shape, seed):
    """
    Return a sequence of uniformly distributed pseudo-random uint32 values.

    The function uses the xorshift pseudo-random number generator, which here
    provides a period length of 2**32 - 1.

    Parameters
    ----------
    shape : int or tuple of int
        Shape of the ndarray to be returned.

    seed : int
        Seed for the random number generator, a positive uint32 value.

    Returns
    -------
    ret : ndarray of uint32
        Array of pseudo-random integers, shaped as specified by `shape`.

    """
    if not _utils.is_scalar_int(seed) or seed < 1 or seed >= 2**32:
        raise ValueError("seed must be a positive integer less than 2**32")

    ret = _np.empty(shape, dtype=_np.double)
    _lib.ordpat_xorshift_rand_double(ret, ret.size, int(seed))

    return ret
