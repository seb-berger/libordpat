/**
 * libordpat - A function library for extracting and encoding
 *             ordinal patterns from time series.
 *
 * For a detailed description of the algorithms used, see the article:
 *
 * [1] Berger S, Kravtsiv A, Schneider G, Jordan D.
 *     Teaching Ordinal Patterns to a Computer.
 *     Entropy. 2019; 21(10):1023.
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

#include "ordpat.h"

#include <stdlib.h>
#include <string.h>

/**
 * function: sequence_length
 *
 * Return the number of ordinal patterns that a time series of length @n_in
 * will yield when transformed using the pattern order @ord and time lag @lag.
 */

static size_t sequence_length(size_t       n_in,
                              unsigned int ord,
                              unsigned int lag)
{
    size_t tmp = (size_t)(ord - 1) * lag;
    return tmp >= n_in ? 0 : n_in - tmp;
}


/**
 * function: check_encoding_args
 *
 * Check whether a combination of encoding parameters @n_in, @n_out, @ord and
 * @lag are valid.
 *
 * Only pattern orders @ord smaller than or equal to @max_ord are accepted.
 *
 * On success, the function stores the resulting number of ordinal patterns
 * in @num_pat, and returns 0. Otherwise, one of the positive error codes
 * defined by @ordpat_error is returned.
 */

static int check_encoding_args(size_t        n_in,
                               size_t        n_out,
                               unsigned int  ord,
                               unsigned int  lag,
                               unsigned int  max_ord,
                               size_t       *num_pat)
{
    if (ord < 2 || ord > max_ord)
        return ORDPAT_ERROR_ORDER_INVALID;

    if (lag < 1)
        return ORDPAT_ERROR_LAG_INVALID;

    *num_pat = sequence_length(n_in, ord, lag);

    if (*num_pat == 0)
        return ORDPAT_ERROR_INPUT_TOO_SHORT;

    if (n_out < *num_pat)
        return ORDPAT_ERROR_OUTPUT_TOO_SHORT;

    return ORDPAT_SUCCESS;
}


/**
 * function: encode_pattern
 *
 * Extract and encode a single ordinal pattern of order @ord from the time
 * series pointed to by @in, using the time lag @tau.
 *
 * The pointer @in must point to an array of at least (@ord - 1) * @lag values.
 * The order @ord must be between 2 and 20. The time lag @lag must be greater
 * than 0.
 *
 * FOR REASONS OF EFFICIENCY, THIS FUNCTION DOES NOT PERFORM ANY ERROR CHECKS!
 *
 * The ordinal pattern is returned as an unsigned integer value between 0 and
 * factorial(@ord) - 1.
 */

static uint64_t encode_pattern(const double *in,
                               unsigned int  ord,
                               unsigned int  lag)
{
    uint64_t code = 0;

    while (--ord) {
        const double *ptr = in;
        unsigned int  pos = ord;

        while (pos--)
            code += *in > *(ptr += lag);

        code *= ord;
        in += lag;
    }

    return code;
}


/**
 * function: add_mp
 *
 * Perform a multi-precision non-negative integer addition.
 *
 * The function adds @op to the value pointed to by @dst. If the addition
 * overflows, a carry value '1' is added to the value pointed to by (@dst + 1).
 * This procedure is repeated until no overflow occurs.
 *
 * Thus, the function interprets @dst as an array of uint64 values, with each
 * array element representing 64 bits of an unsigned integer of arbitrary size.
 * Independent of the endianness of the underlying hardware, the function uses
 *
 *     @dst[0] for the bits   0 to  63,
 *     @dst[1] for the bits  64 to 127,
 *     @dst[2] for the bits 128 to 191
 *
 * and so forth.
 *
 * It is the caller's responsibility that the array pointed to by @dst is big
 * enough to store the overall result. Unbounded memory access will otherwise
 * occur.
 */

static void add_mp(uint64_t *dst,
                   uint64_t  op)
{
    while (op) {
        uint64_t tmp = *dst;
        *dst += op;
        op = *dst++ < tmp;
    }
}


/**
 * function: subtract_mp
 *
 * Perform a multi-precision non-negative integer subtraction.
 *
 * The function subtracts @op from the value pointed to by @dst. If the
 * subtraction underflows, a carry value '1' is subtracted from the value
 * pointed to by (@dst + 1). This procedure is repeated until no underflow
 * occurs.
 *
 * Thus, the function interprets @dst as an array of uint64 values, with each
 * array element representing 64 bits of an unsigned integer of arbitrary size.
 * Independent of the endianness of the underlying hardware, the function uses
 *
 *     @dst[0] for the bits   0 to  63,
 *     @dst[1] for the bits  64 to 127,
 *     @dst[2] for the bits 128 to 191
 *
 * and so forth.
 *
 * It is the caller's responsibility that the result of the operation is
 * non-negative. Unbounded memory access will otherwise occur.
 */

static void subtract_mp(uint64_t *dst,
                        uint64_t  op)
{
    while (op) {
        uint64_t tmp = *dst;
        *dst -= op;
        op = *dst++ > tmp;
    }
}


/**
 * function: multiply_mp
 *
 * Perform a multi-precision non-negative integer multiplication.
 *
 * The function multiplies @op with the (@n times 64) bit wide integer pointed
 * to by @dst, and stores the result in @dst.
 *
 * Thus, the function interprets @dst as an array of uint64 values, with each
 * array element representing 64 bits of an unsigned integer of arbitrary size.
 * Independent of the endianness of the underlying hardware, the function uses
 *
 *     @dst[0] for the bits   0 to  63,
 *     @dst[1] for the bits  64 to 127,
 *     @dst[2] for the bits 128 to 191
 *
 * and so forth.
 *
 * The result will silently be truncated to its lowest (@n times 64) bit, that
 * is, overflow will wrap around in a 'modulo-like' manner.
 *
 */

static void multiply_mp(uint64_t *dst,
                        size_t    n,
                        uint32_t  op)
{
    dst  += n - 1;
    *dst *= op;

    while (--n) {
        uint64_t high;

        --dst;

        high = (*dst >> 32) * op;
        *dst = (*dst & 0xFFFFFFFF) * op;

        add_mp(dst,     high << 32);
        add_mp(dst + 1, high >> 32);
    }
}


/**
 * function: swap
 *
 * Swap the value pointed to by @x with the value pointed to by @y.
 */

static void swap(double *x,
                 double *y)
{
    double tmp = *x;
    *x = *y;
    *y = tmp;
}


/**
 * function: reverse
 *
 * Reverse the order of the first @len elements in the array pointed to
 * by @seq.
 */

static void reverse(double *seq,
                    size_t  len)
{
    double *beg = seq;
    double *end = beg + len - 1;

    len /= 2;

    while (len--)
        swap(beg++, end--);
}


/**
 * function: next_perm
 *
 * Permute the first @len elements of the array pointed to by @tuple, such that
 * the result is the 'next' permutation with regard to the lexicographic order
 * of all the factorial(@len) permutations of @len elements.
 *
 * If no such permutation exists, the function does not alter the tuple. For
 * example, there are 6 tuples of 3 elements, and the function turns
 *
 *     (1, 2, 3) into (1, 3, 2),
 *     (1, 3, 2) into (2, 1, 3),
 *     (2, 1, 3) into (2, 3, 1),
 *     (2, 3, 1) into (3, 1, 2),
 *     (3, 1, 2) into (3, 2, 1),
 *
 * whereas (3, 2, 1) remains (3, 2, 1).
 *
 * The function returns 0 if a new permutation was created, or 1 otherwise.
 * The rationale is that one can obtain all permutations of a given length by
 * starting with a tuple sorted in ascending order, and repeatedly calling
 * @next_perm until 1 is returned.
 *
 * The function requires that all elements in the array @tuple are pairwise
 * distinct. Unexpected behaviour may otherwise occur.
 */

static int next_perm(double *tuple,
                     size_t  len)
{
    size_t pos1, pos2;

    pos1 = len;
    while (--pos1 && tuple[pos1] < tuple[pos1 - 1]);

    if (pos1 == 0)
        return 1;

    --pos1;

    pos2 = len;
    while (--pos2 != pos1 && tuple[pos2] < tuple[pos1]);

    swap(tuple + pos1, tuple + pos2);
    reverse(tuple + pos1 + 1, len - pos1 - 1);
    return 0;
}


/**
 * function: xorshift32
 *
 * Return the next xorshifted pseudo-random number following @x, using the
 * hard-coded parameters below.
 */

static uint32_t xorshift32(uint32_t x)
{
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    return x;
}


/****** THE EXPORTED LIBRARY FUNCTIONS BELOW ARE DOCUMENTED IN ORDPAT.H ******/

ORDPAT_API uint64_t ordpat_factorial(unsigned int x)
{
    /* factorials from 0! to 20! */
    static const uint64_t table[21] = {
        1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800,
        479001600, 6227020800, 87178291200, 1307674368000, 20922789888000,
        355687428096000, 6402373705728000, 121645100408832000,
        2432902008176640000
    };

    return x <= 20 ? table[x] : 0;
}


ORDPAT_API size_t ordpat_pattern_width(unsigned int ord)
{
    /**
     * In essence, the function calculates factorial(@ord) - 1, and
     * then determines the width of the result in multiples of 64 bit.
     */

    uint64_t *fact;  /* storage space for multi-precision factorial */
    size_t    alloc;

    if (ord < 2 || ord > 255) /* Orders > 255 are not supported. */
        return 0;

    alloc   = 2;
    fact    = calloc(alloc, sizeof(uint64_t));
    fact[0] = ord;

    while (--ord > 1) {
        multiply_mp(fact, alloc, ord);

        /* Always keep an extra 64 bit of zeros to detect overflows. */
        if (fact[alloc - 1] != 0) {
            ++alloc;
            fact = realloc(fact, alloc * sizeof(uint64_t));
            fact[alloc - 1] = 0;
        }
    }

    subtract_mp(fact, 1); /* equals: factorial(@ord) - 1 */

    /* Count the number of 64-bit words used. */
    while (fact[alloc - 1] == 0)
        --alloc;

    free(fact);
    return alloc;
}


ORDPAT_API int ordpat_create_lookup_table(uint64_t     *table,
                                          size_t        len,
                                          unsigned int  ord)
{
    double    tuple[10] = {2, 4, 6, 8, 10, 12, 14, 16, 18, 20};
    uint64_t *ptr       = table;
    size_t    fact;
    size_t    num_tails;

    if (ord < 2 || ord > 10)
        return ORDPAT_ERROR_ORDER_INVALID;

    num_tails = ordpat_factorial(ord - 1);
    fact = num_tails * ord;

    if (fact * ord > len)
        return ORDPAT_ERROR_INPUT_TOO_SHORT;

    /* Generate and encode all possible pattern combinations. */
    while (num_tails--) {
        size_t num_ranks = ord;

        while (num_ranks--) {
            tuple[ord - 1] = 2 * num_ranks + 1;
            *ptr++ = encode_pattern(tuple, ord, 1);
        }

        next_perm(tuple, ord - 1);
    }

    while (--ord) {
        memcpy(ptr, table, fact * sizeof(uint64_t));
        ptr += fact;
    }

    return ORDPAT_SUCCESS;
}


ORDPAT_API int ordpat_encode_plain(const double *in,
                                   size_t        n_in,
                                   uint64_t     *out,
                                   size_t        n_out,
                                   unsigned int  ord,
                                   unsigned int  lag)
{
    size_t len;
    int    ret;

    ret = check_encoding_args(n_in, n_out, ord, lag, 20, &len);

    if (ret != ORDPAT_SUCCESS)
        return ret;

    /* Plain and simple, hence the name... */
    while (len--)
        *out++ = encode_pattern(in++, ord, lag);

    return ORDPAT_SUCCESS;
}


ORDPAT_API int ordpat_encode_overlap(const double *in,
                                     size_t        n_in,
                                     uint64_t     *out,
                                     size_t        n_out,
                                     unsigned int  ord,
                                     unsigned int  lag)
{
    const double *next;
    uint8_t       stack_mem[20 * 25] = {0};
    uint8_t      *mem;
    uint8_t      *buf;
    size_t        len;
    size_t        pos;
    int           ret;

    ret = check_encoding_args(n_in, n_out, ord, lag, 20, &len);

    if (ret != ORDPAT_SUCCESS)
        return ret;

    /**
     * The size of @stack_mem should be enough for most applications.
     * Otherwise, we prepare some memory on the heap.
     */
    if (lag * ord <= sizeof(stack_mem))
        mem = stack_mem;
    else
        mem = calloc(lag * ord, 1);

    --ord; /* makes the code a little simpler */

    /* Extract inversion counts for the first @lag ordinal patterns. */
    buf = mem;
    for (pos = 0; pos < lag; ++pos) {
        unsigned int i, j;

        for (i = 0; i < ord - 1; ++i) {
            for (j = i + 1; j < ord; ++j) {
                *buf += (in[pos + i * lag] >
                         in[pos + j * lag]);
            }
            ++buf;
        }
    }

    buf  = mem;
    pos  = 0;
    next = in + ord * lag;

    while (len--) {
        unsigned int i;

        /* No need to store the leftmost inversion count in the buffer. */
        *out = ord * (*buf + (*in > *next));

        /**
         * This is the 'vanilla' encoding approach. The line above, and the
           two lines below the for loop are simplifications for the leftmost
           and rightmost tuple entries.
         */
        for (i = 1; i < ord - 1; ++i) {
            /* Shift buffer and update inversion count. */
            *buf = *(buf + 1) + (in[i * lag] > *next);

            /* Accumulate inversion count onto output code. */
            *out += *buf++;
            *out *= ord - i;
        }

        /* The rightmost inversion count does not use previous results. */
        *buf = in[i * lag] > *next;

        /* Also, no multiplication (by 1) is required. */
        *out += *buf;

        /* Step ahead. */
        ++in;
        ++out;
        ++buf;
        ++next;

        /* Rewind circular buffer if necessary. */
        if (++pos == lag) {
            buf = mem;
            pos = 0;
        }
    }

    /* Release heap memory if used. */
    if (mem != stack_mem)
        free(mem);

    return ORDPAT_SUCCESS;
}


ORDPAT_API int ordpat_encode_overlap_mp(const double *in,
                                        size_t        n_in,
                                        uint64_t     *out,
                                        size_t        n_out,
                                        unsigned int  ord,
                                        unsigned int  lag)
{
    const double *next;
    uint8_t       stack_mem[255 * 25] = {0};
    uint8_t      *mem;
    uint8_t      *buf;
    size_t        len;
    size_t        width;
    size_t        pos;
    int           ret;

    ret = check_encoding_args(n_in, n_out, ord, lag, 255, &len);

    if (ret != ORDPAT_SUCCESS)
        return ret;

    width = ordpat_pattern_width(ord);

    if (n_out < len * width)
        return ORDPAT_ERROR_OUTPUT_TOO_SHORT;

    /**
     * The size of @stack_mem should be enough for most applications.
     * Otherwise, we prepare some memory on the heap.
     */
    if (lag * ord > sizeof(stack_mem))
        mem = calloc(lag * ord, 1);
    else
        mem = stack_mem;

    --ord; /* makes the code a little simpler */

    /* Extract inversion counts for the first @lag ordinal patterns. */
    buf = mem;
    for (pos = 0; pos < lag; ++pos) {
        unsigned int i, j;

        for (i = 0; i < ord - 1; ++i) {
            for (j = i + 1; j < ord; ++j) {
                *buf += (in[pos + i * lag] >
                         in[pos + j * lag]);
            }
            ++buf;
        }
    }

    buf   = mem;
    pos   = 0;
    next  = in + ord * lag;

    /**
     * This while loop works exactly like the one in @encode_ordpat_overlap,
     * but uses multiple-precision arithmetic to support pattern orders > 20.
     */
    while (len--) {
        unsigned int i;

        /* No need to store the leftmost inversion count in the buffer. */
        *out = *buf + (*in > *next);
        multiply_mp(out, width, ord);

        for (i = 1; i < ord - 1; ++i) {
            /* Shift buffer and update inversion count. */
            *buf = *(buf + 1) + (in[i * lag] > *next);

            /* Accumulate inversion count onto output code. */
            add_mp(out, *buf++);
            multiply_mp(out, width, ord - i);
        }

        /* The rightmost inversion count does not use previous results. */
        *buf = in[i * lag] > *next;

        /* Also, no multiplication (by 1) is required. */
        add_mp(out, *buf);

        /* Step ahead. */
        ++in;
        ++buf;
        ++next;
        out += width;

        /* Rewind circular buffer if necessary. */
        if (++pos == lag) {
            pos = 0;
            buf = mem;
        }
    }

    /* Release heap memory if used. */
    if (mem != stack_mem)
        free(mem);

    return ORDPAT_SUCCESS;
}


ORDPAT_API int ordpat_encode_lookup(const double *in,
                                    size_t        n_in,
                                    uint64_t     *out,
                                    size_t        n_out,
                                    unsigned int  ord,
                                    unsigned int  lag,
                                    uint64_t     *tab)
{
    size_t len;
    size_t pos;
    int    ret;

    ret = check_encoding_args(n_in, n_out, ord, lag, 10, &len);

    if (ret != ORDPAT_SUCCESS)
        return ret;

    /**
     * Obtain the first @lag ordinal patterns directly.
     *
     * (This is the reason why @tab MUST point to a lookup table created by
     *  the function @ordpat_create_lookup_table, and no other bijection
     *  between ordinal patterns and numerical codes is supported.)
     */
    for (pos = 0; pos < lag; ++pos)
        *out++ = encode_pattern(in++, ord, lag);

    /* The remaining patterns can be determined using the lookup table: */
    for (pos = lag; pos < len; ++pos) {
        unsigned int count = 0;
        unsigned int i;

        /* Get the rightmost inversion count. */
        for (i = 0; i < ord - 1; ++i)
            count += in[i * lag] > in[(ord - 1) * lag];

        /**
         * Look up the result, which depends on the inversion count, as well
         * as the ordinal pattern @lag time steps earlier.
         */
        *out = tab[*(out - lag) * ord + count];

        /* Step ahead. */
        ++in;
        ++out;
    }

    return ORDPAT_SUCCESS;
}


ORDPAT_API void ordpat_xorshift_rand_uint32(uint32_t *dst,
                                            size_t    len,
                                            uint32_t  seed)
{
    while (len--)
        *dst++ = seed = xorshift32(seed);
}


ORDPAT_API void ordpat_xorshift_rand_double(double  *dst,
                                            size_t   len,
                                            uint32_t seed)
{
    while (len--)
        *dst++ = seed = xorshift32(seed);
}
