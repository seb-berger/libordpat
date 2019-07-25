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
libordpat - A library for extracting ordinal patterns from time series.

Example of usage
----------------
>>> import numpy as np
>>> import ordpat
>>> x = np.random.randn((4, 10000))
>>> y = ordpat.ordpat(x, 5, 2)
>>> y.shape
(4, 9992)

For a detailed description of the algorithms used, see the article:

[1] Berger S, Kravtsiv A, Schneider G, Jordan D.
    Teaching Ordinal Patterns to a Computer.
    Entropy. 2019; xx(xx):xxx.
    https://doi.org/xx.xxxx/xxxxxxxxx

"""
from __future__ import division as _division
import numpy as _np

import encoders as native
import utils as _utils

# Try to import libordpat_wrapper, which will fail if the libordpat shared
# library cannot be dynamically loaded. This exception is handled gracefully,
# such that the Python module is still useful.
try:
    import library_wrapper as library
    _library_loaded = True
except OSError:
    # Create a dummy module
    class library:
        """
        WARNING: The libordpat C library was not found.
        """
        encode_lookup     = None
        encode_overlap    = None
        encode_overlap_mp = None
        encode_plain      = None

    _library_loaded = False
    print("WARNING: Could not find the C library.")

# This is used to keep the last lookup table in memory, such that multiple
# calls using the same pattern order do not have to recreate a lookup table
# each time.
_lookup_buf = {'order': None, 'table': None}


# This data structure describes all the algorithms available, including their
# maximum pattern order, respective function call, as well as a flag that
# indicates whether the libordpat C library is required.
_algorithms = {
    'lookup' : {
        'max_order' : 10,
        'c_lib'     : False,
        'fcn'       : lambda x, ord, lag, trustme : \
                      native.encode_lookup(x, ord, lag, _lookup_buf['table'],
                                           trustme)
        },

    'lookup_c' : {
        'max_order' : 10,
        'c_lib'     : True,
        'fcn'       : lambda x, ord, lag, trustme: \
                      library.encode_lookup(x, ord, lag, _lookup_buf['table'],
                                            trustme)
        },

    'overlap' : {
        'max_order' : 20,
        'c_lib'     : False,
        'fcn'       : native.encode_overlap
        },

    'overlap_c' : {
        'max_order' : 20,
        'c_lib'     : True,
        'fcn'       : library.encode_overlap
        },

    'overlap_mp_c' : {
        'max_order' : 255,
        'c_lib'     : True,
        'fcn'       : library.encode_overlap_mp
        },

    'plain' : {
        'max_order' : 20,
        'c_lib'     : False,
        'fcn'       : native.encode_plain
        },

    'plain_c' : {
        'max_order' : 20,
        'c_lib'     : True,
        'fcn'       : library.encode_plain
        },

    'vectorised' : {
        'max_order' : 20,
        'c_lib'     : False,
        'fcn'       : native.encode_vectorised
        },
    }


def _update_lookup_table(ord):
    if _lookup_buf['order'] != ord:
        if _library_loaded:
            _lookup_buf['table'] = library.create_lookup_table(ord);
        else:
            _lookup_buf['table'] = native.create_lookup_table(ord);

        _lookup_buf['order'] = ord;


def _check_algorithm(ord, alg):
    if alg is None:
        if _library_loaded:
            if ord <= 20:
                alg = "overlap_c"
            else:
                alg = "overlap_mp_c"
        elif ord <= 20:
            alg = "vectorised"
        else:
            raise ValueError("pattern orders > 20 require the C library")

    elif alg not in _algorithms:
        if isinstance(alg, str):
            raise ValueError("unknown algorithm: '{}'".format(alg))
        else:
            raise TypeError("algorithm must be a string")

    elif ord > _algorithms[alg]['max_order']:
        raise ValueError("{} algorithm does not support orders > {}".format(
                         alg, _algorithms[alg]['max_order']))

    elif _algorithms[alg]['c_lib'] and not _library_loaded:
        raise ValueError("{} algorithm requires the " \
                         "C library".format(alg))

    return alg


def ordpat(x, ord, lag=1, algorithm=None, axis=-1):
    """
    Return a sequence of ordinal patterns in numerical representation.

    The function extracts the consecutive ordinal patterns of order `ord`
    from the (possibly multi-variate) time series contained in `x`, using the
    time lag `lag` and the encoding algorithm `algorithm`.

    The ordinal patterns are then encoded into numerical values ranging
    from 0 to factorial(ord) - 1.

    Parameters
    ----------
    x : array_like
        Time series data to be encoded. Can be multi-dimensional, and must be
        convertible to an ndarray of data type 'double'.

    ord : int
        Order of the encoding, an integer greater than 2. The upper bound is
        algorithm-dependent (see below).

    lag : int, optional
        Time lag used for creating embedding vectors. Must be a positive
        integer.

    algorithm : {'plain', 'plain_c', 'overlap', 'overlap_c', 'overlap_mp_c',
                 'lookup', 'lookup_c', 'vectorised'}, optional
        One of the encoding algorithms described below. By default, the
        algorithm is selected automatically.

    axis : int, optional
        The axis of the input data along which the extraction shall be
        performed. By default, the last axis is used.

    Returns
    -------
    pat : ndarray of uint64
        Ordinal patterns in numerical representation and stored as integers
        in uint64 format. The shape of the `pat` array is equal to the input
        array `x`, except for the axis on which the encoding is performed.
        For this axis, it holds that

          ``pat.shape[axis] == x.shape[axis] - (ord - 1) * lag``.

        The 'overlap_mp_c' algorithm represents ordinal patterns by vectors
        of uint64 values. Thus, an extra dimension is appended to `pat` in
        this case. The size of this dimension depends on the pattern order.

    Details
    -------
    The `ordpat` function provides various encoding algorithms, all yielding
    the same numerical result. The algorithms differ in run-time efficiency,
    as well as in the maximum pattern order supported.

    Possible values for `algorithm` are:

    --------------------------------------------------------------------------
    ALGORITHM    MAX ORD    NEEDS C LIBRARY    REMARKS
    --------------------------------------------------------------------------
    lookup            10                 no
    overlap           20                 no
    plain             20                 no    Very slow.
    vectorised        20                 no    Best native choice.
    --------------------------------------------------------------------------
    lookup_c          10                yes    Fast for small orders.
    overlap_c         20                yes    Good general-purpose algorithm.
    overlap_mp_c     255                yes    Supports orders > 20.
    plain_c           20                yes
    --------------------------------------------------------------------------

    The 'overlap_mp_c' algorithm supports very high pattern orders.
    When using this algorithm, the last dimension of the returned array
    indexes the 64-bit wide elements of each single ordinal pattern.
    For example,

    >>>  ordpat(np.zeros((8, 1000)), 21, 1, 'overlap_mp_c')

    returns a uint64-valued array of shape (8, 980, 2), because a sequence
    of length 1000 has 980 ordinal patterns of order 21, whereas storing
    each such pattern requires two 64-bit values.

    """
    x = _utils.prepare_input_data(x, ord, lag, axis)
    alg = _check_algorithm(ord, algorithm)

    enc_fnc = _algorithms[alg]['fcn']

    if alg in ("lookup", "lookup_c"):
        _update_lookup_table(ord)

    if alg == "vectorised":
        return enc_fnc(x, ord, lag, axis, True)
    else:
        return _np.apply_along_axis(enc_fnc, axis, x, ord, lag, True)
