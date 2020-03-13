/**
 * libordpat - A function library for extracting and encoding
 *             ordinal patterns from time series.
 *
 * For a detailed description of the algorithms used, see the article:
 *
 * [1] Berger S, Kravtsiv A, Schneider G, Jordan D.
 *     Teaching Ordinal Patterns to a Computer: Efficient Encoding
 *     Algorithms Based on the Lehmer Code. Entropy. 2019; 21(10):1023.
 *     https://doi.org/10.3390/e21101023
 */

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

#ifndef ORDPAT_H
#define ORDPAT_H

#if defined(_WIN32)
    #if defined(BUILD_LIBORDPAT)
        #if defined(__GNUC__)
            #define ORDPAT_API __attribute__ ((dllexport))
        #else
            #define ORDPAT_API __declspec(dllexport)
        #endif
    #elif defined(BUILD_LIBORDPAT_STATIC)
            #define ORDPAT_API
    #else
        #if defined(__GNUC__)
            #define ORDPAT_API __attribute__ ((dllimport))
        #else
            #define ORDPAT_API __declspec(dllimport)
        #endif
    #endif
#elif defined(__unix) || defined(__unix__) || defined(__APPLE_CC__)
    #if defined(BUILD_LIBORDPAT)
        #define ORDPAT_API __attribute__ ((visibility ("default")))
    #else
        #define ORDPAT_API
    #endif
#else
    #error Unknown platform.
#endif

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>
#include <stdint.h>


/**
 * enumeration: ordpat_error
 *
 * Error codes returned by the encoding functions.
 */

enum ordpat_error {
    ORDPAT_SUCCESS                = 0,
    ORDPAT_ERROR_ORDER_INVALID    = 1,
    ORDPAT_ERROR_LAG_INVALID      = 2,
    ORDPAT_ERROR_INPUT_TOO_SHORT  = 3,
    ORDPAT_ERROR_OUTPUT_TOO_SHORT = 4
};


/**
 * function: ordpat_factorial
 *
 * Return the factorial of @x for @x < 21, or 0 otherwise.
 */

ORDPAT_API uint64_t ordpat_factorial(unsigned int x);


/**
 * function: ordpat_pattern_width
 *
 * Return the number of 64-bit words required to store an ordinal pattern
 * of order @ord, where 1 < @ord < 255, or 0 otherwise.
 *
 * The function assumes that each pattern of order @ord shall be represented
 * by a unique non-negative integer in {0, 1, ..., factorial(@ord) - 1}.
 */

ORDPAT_API size_t ordpat_pattern_width(unsigned int ord);


/**
 * function: ordpat_create_lookup_table
 *
 * Create a lookup table for the 'lookup' encoding algorithm (see [1]).
 *
 * Orders 1 < @ord < 11, and any time lag @lag > 0 are supported.
 *
 * The array @table needs to be pre-allocated to @len = factorial(@ord) * @ord
 * elements. On return of the function, the @table will contain values between
 * 0 and factorial(@ord) - 1 exclusively. In particular, it is initialised
 * such that
 *
 *     @n_1 = @table[@n_0 * @ord + @inv_count],
 *
 * where @n_0 is the ordinal pattern of a tuple (x_0, x_1, ..., x_{@ord - 1}),
 * @n_1 is the ordinal pattern of the subsequent tuple (x_1, x_2, ..., x_@ord),
 * and @inv_count is the number of left inversions of x_@ord with regard to
 * x_1, x_2, ..., x_{@ord - 1}.
 *
 * The function returns 0 on success, or one of the positive error codes
 * defined by @ordpat_error.
 */

ORDPAT_API int ordpat_create_lookup_table(uint64_t     *table,
                                          size_t        len,
                                          unsigned int  ord);


/**
 * function: ordpat_encode_plain
 *
 * Transform the double-valued time series @in into a sequence of ordinal
 * patterns using the 'plain' algorithm (see [1]).
 *
 * Orders 1 < @ord < 21, and any time lag @lag > 0 are supported.
 *
 * Ordinal patterns of order @ord are extracted from the time series @in of
 * length @n_in, using the time lag @lag. For a given @ord, there are
 * factorial(@ord) different ordinal patterns. A non-negative integer in
 * {0, 1, ..., factorial(@ord) - 1} is used to uniquely represent each such
 * pattern.
 *
 * The resulting numerical codes are successively stored in @out, which must
 * be pre-allocated to an appropriate size. The number of elements available
 * in @out has to be provided in @n_out, and must be at least @n_out =
 * @n_in - (@ord - 1) * @lag.
 *
 * The function returns 0 on success, or one of the positive error codes
 * defined by @ordpat_error.
 */

ORDPAT_API int ordpat_encode_plain(const double *in,
                                   size_t        n_in,
                                   uint64_t     *out,
                                   size_t        n_out,
                                   unsigned int  ord,
                                   unsigned int  lag);


/**
 * function: ordpat_encode_overlap
 *
 * Transform the double-valued time series @in into a sequence of ordinal
 * patterns using the 'overlap' algorithm (see [1]).
 *
 * Orders 1 < @ord < 21, and any time lag @lag > 0 are supported.
 *
 * Ordinal patterns of order @ord are extracted from the time series @in of
 * length @n_in, using the time lag @lag. For a given @ord, there are
 * factorial(@ord) different ordinal patterns. A non-negative integer in
 * {0, 1, ..., factorial(@ord) - 1} is used to uniquely represent each such
 * pattern.
 *
 * The resulting numerical codes are successively stored in @out, which must
 * be pre-allocated to an appropriate size. The number of elements available
 * in @out has to be provided in @n_out, and must be at least @n_out =
 * @n_in - (@ord - 1) * @lag.
 *
 * The function returns 0 on success, or one of the positive error codes
 * defined by @ordpat_error.
 */

ORDPAT_API int ordpat_encode_overlap(const double *in,
                                     size_t        n_in,
                                     uint64_t     *out,
                                     size_t        n_out,
                                     unsigned int  ord,
                                     unsigned int  lag);


/**
 * function: ordpat_encode_overlap_mp
 *
 * Transform the double-valued time series @in into a sequence of ordinal
 * patterns using the 'overlap_mp' algorithm (see [1]).
 *
 * Orders 1 < @ord < 255, and any time lag @lag > 0 are supported.
 *
 * Ordinal patterns of order @ord are extracted from the time series @in of
 * length @n_in, using the time lag @lag. For a given @ord, there are
 * factorial(@ord) different ordinal patterns. A non-negative integer in
 * {0, 1, ..., factorial(@ord) - 1} is used to uniquely represent each such
 * pattern.
 *
 * The resulting numerical codes are successively stored in @out, which must
 * be pre-allocated to an appropriate size. For @ord > 20, the width of a
 * single ordinal pattern may exceed 64 bit. Therefore, each ordinal pattern
 * is represented by a tuple of consecutive uint64 values.
 *
 * The number of elements available in @out has to be provided in @n_out,
 * and must be at least @n_out = @width * (@n_in - (@ord - 1) * @lag), where
 * @width is the number of 64-bit values that make up one ordinal pattern.
 * The value of @width for a given @ord can be obtained by calling the
 * function @ordpat_pattern_width.
 *
 * The function returns 0 on success, or one of the positive error codes
 * defined by @ordpat_error.
 */

ORDPAT_API int ordpat_encode_overlap_mp(const double *in,
                                        size_t        n_in,
                                        uint64_t     *out,
                                        size_t        n_out,
                                        unsigned int  ord,
                                        unsigned int  lag);


/**
 * function: ordpat_encode_lookup
 *
 * Transform the double-valued time series @in into a sequence of ordinal
 * patterns using the 'lookup' algorithm (see [1]).
 *
 * Orders between 1 < @ord < 11, and any time lag @lag > 0 are supported.
 *
 * Ordinal patterns of order @ord are extracted from the time series @in of
 * length @n_in, using the time lag @lag. For a given @ord, there are
 * factorial(@ord) different ordinal patterns. A non-negative integer in
 * {0, 1, ..., factorial(@ord) - 1} is used to uniquely represent each such
 * pattern.
 *
 * The resulting numerical codes are successively stored in @out, which must
 * be pre-allocated to an appropriate size. The number of elements available
 * in @out has to be provided in @n_out, and must be at least @n_out =
 * @n_in - (@ord - 1) * @lag.
 *
 * The 'lookup' algorithms relies on a lookup table @tab, which has to be
 * created by calling @ordpat_create_lookup_table. The idea is to create
 * the lookup table only once, then reuse it for multiple invocations of the
 * encoding function (instead of keeping the table internal to the encoding
 * function, and recreating it each time the function is called).
 *
 * WARNING: Note that @ordpat_encode_lookup does NOT support user-defined
 * lookup tables. Any other than the table returned by @ordpat_create_lookup
 * will result in erroneous results (at best).
 *
 * The function returns 0 on success, or one of the positive error codes
 * defined by @ordpat_error.
 */

ORDPAT_API int ordpat_encode_lookup(const double *in,
                                    size_t        n_in,
                                    uint64_t     *out,
                                    size_t        n_out,
                                    unsigned int  ord,
                                    unsigned int  lag,
                                    uint64_t     *tab);


/**
 * function: ordpat_xorshift_rand_uint32
 *
 * Create @len pseudo-random numbers between 0 and 2^32 - 1, and store them
 * in the array @dst, which needs to be sized appropriately.
 *
 * Numbers are produced using the xorshift random number generator, which is
 * initialised to @seed before the first value is obtained. Thus, the first
 * element in @dst is the number immediately following @seed.
 */

ORDPAT_API void ordpat_xorshift_rand_uint32(uint32_t *dst,
                                            size_t    len,
                                            uint32_t  seed);


/**
 * function: ordpat_xorshift_rand_double
 *
 * A convenience function that does the same as @xorshift_rand_uint32, but
 * stores its result in an array of double-precision floating point values
 * instead.
 */

ORDPAT_API void ordpat_xorshift_rand_double(double  *dst,
                                            size_t   len,
                                            uint32_t seed);

#ifdef __cplusplus
}
#endif

#endif /* ORDPAT_H */
