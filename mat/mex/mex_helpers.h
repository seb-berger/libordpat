/**
 * Copyright (c) 2019, Sebastian Berger.
 *
 * Klinikum rechts der Isar der
 * Technischen Universitaet Muenchen
 * Munich, Germany
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in
 *       the documentation and/or other materials provided with the
 *       distribution.
 *     * Neither the names of the copyright holders nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR OR
 * THE KLINIKUM RECHTS DER ISAR BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
 * THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
 * DAMAGE.
 */

#ifndef MEX_HELPERS_H
#define MEX_HELPERS_H

#include <stdio.h>
#include <stdint.h>
#include <mex.h>

/**
 * function: is_mw_size
 *
 * Return true if @d can safely be converted to 'mwSize', false otherwise.
 */

static inline bool is_mw_size(double d)
{
    /**
     * When compiled using the parameter '--compatibleArrayDims' on MATLAB,
     * mwSize is a typedef for signed int. We want the function to
     * consistently return false for negative values, though.
     */

    return mxIsFinite(d) && d >= 0 && (double)(mwSize)(d) == d;
}


/**
 * function: is_integer_64
 *
 * Return true if @d can safely be converted to 'int64_t', false otherwise.
 */

static inline bool is_integer_64(double d)
{
    return mxIsFinite(d) && (double)(int64_t)(d) == d;
}


/**
 * function: is_unsigned_integer_32
 *
 * Return true if @d can safely be converted to 'uint32_t', false otherwise.
 */

static inline bool is_unsigned_integer_32(double d)
{
    return mxIsFinite(d) && (double)(uint32_t)(d) == d;
}


/**
 * function: narginchk
 *
 * Raise an error if @n_in is smaller than @min or greater than @max.
 * Otherwise return to caller.
 */

static inline void narginchk(int n_in,
                             int min,
                             int max)
{
    if (n_in < min)
        mexErrMsgTxt("Not enough input arguments.");

    if (n_in > max)
        mexErrMsgTxt("Too many input arguments.");
}


/**
 * function: nargoutchk
 *
 * Raise an error if @n_in is smaller than @min or greater than @max.
 * Otherwise return to caller.
 */

static inline void nargoutchk(int n_out,
                              int min,
                              int max)
{
    if (n_out < min)
        mexErrMsgTxt("Not enough output arguments.");

    if (n_out > max)
        mexErrMsgTxt("Too many output arguments.");
}


/**
 * function: mx_to_double_array
 *
 * Return a const pointer to the double-precision data in the mxArray @in.
 * Raise an error if @in does not contain real-valued, double-precision data.
 */

static inline const double* mx_to_double_array(const mxArray *in)
{
    if (!mxIsDouble(in) || mxIsComplex(in))
        mexErrMsgTxt("X must be in double-precision float format.");

    return mxGetPr(in);
}


/**
 * function: mx_to_pattern_order
 *
 * Return the unsigned integer pattern order contained in the mxArray @in.
 * Raise an error if @in is not a valid pattern order, that is, does not
 * contain an integral number between 2 and 255.
 */

static inline unsigned int mx_to_pattern_order(const mxArray *in)
{
    const char *msg = "ORD must be an integer between 2 and 255.";

    double ord;
    bool valid;

    valid = mxIsNumeric(in) && !mxIsComplex(in) &&
            mxGetNumberOfElements(in) == 1;

    if (!valid)
        mexErrMsgTxt(msg);

    ord = mxGetScalar(in);

    if (!is_unsigned_integer_32(ord) || ord < 2 || ord > 255)
        mexErrMsgTxt(msg);

    return ord;
}


/**
 * function: mx_to_time_lag
 *
 * Return the unsigned integer time lag contained in the mxArray @in.
 * Raise an error if @in is not a valid time lag, that is, does not
 * contain an integral number between 1 and 2^32 - 1.
 */

static inline unsigned int mx_to_time_lag(const mxArray *in)
{
    const char *msg = "LAG must be an integer between 1 and 2^32 - 1.";

    double lag;
    bool valid;

    valid = mxIsNumeric(in) && !mxIsComplex(in) &&
            mxGetNumberOfElements(in) == 1;

    if (!valid)
        mexErrMsgTxt(msg);

    lag = mxGetScalar(in);

    if (!is_unsigned_integer_32(lag) || lag < 1)
        mexErrMsgTxt(msg);

    return lag;
}


/**
 * function: sequence_length
 *
 * Return the number of ordinal patterns that would result from encoding
 * a time series of length @n_input, using a pattern order @ord and a time
 * lag @lag.
 *
 * Raise an error if the time series is too short to hold any ordinal patterns
 * of order @ord and time lag @tau.
 *
 * Note that @ord and @lag must be VALID encoding parameters. Otherwise, the
 * behaviour of the function is undefined.
 */

static inline size_t sequence_length(size_t       n_input,
                                     unsigned int ord,
                                     unsigned int lag)
{
    size_t tmp = ord - 1;
    tmp *= lag;

    if (tmp >= n_input)
        mexErrMsgTxt("Input data X too short for selected ORD and LAG.");

    return n_input - tmp;
}


/**
 * function:catch_libordpat_error
 *
 * Raise an error if @ret is not zero. This is meant to handle return codes
 * of library calls consistently in one place.
 */

static inline void catch_libordpat_error(int ret)
{
    char msg[255] = {0};

    if (ret == 0)
        return;

    snprintf(msg, sizeof(msg), "Encoding function failed with "
             "error code %d. This is a bug!", ret);

    mexErrMsgTxt(msg);
}

#endif /* MEX_HELPERS_H */
