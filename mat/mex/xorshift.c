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

#include <mex.h>
#include <stdint.h>

#include "ordpat.h"
#include "mex_helpers.h"


static mxArray* allocate_output(const mxArray *arg)
{
    mxArray      *out    = NULL;
    mwSize       *sizes  = NULL;
    const double *p_arg  = NULL;
    size_t        n_dims = 0;
    size_t        i      = 0;

    if (mxIsNumeric(arg) == 0
    ||  mxIsComplex(arg) == 1
    ||  mxIsDouble(arg)  == 0
    ||  mxGetNumberOfDimensions(arg) != 2
    ||  mxGetDimensions(arg)[0] != 1)
        mexErrMsgTxt("Size vector must be a row vector of real elements.");

    p_arg  = mxGetPr(arg);
    n_dims = mxGetDimensions(arg)[1];

    sizes = mxCalloc(n_dims, sizeof(mwSize));

    for (i = 0; i < n_dims; ++i) {
        double d = p_arg[i] < 0 ? 0 : p_arg[i];

        if (is_mw_size(d) == false)
            mexErrMsgTxt("Invalid dimensions.");

        sizes[i] = d;
    }

    out = mxCreateNumericArray(n_dims, sizes, mxDOUBLE_CLASS, mxREAL);
    mxFree(sizes);

    return out;
}


static uint32_t read_seed(const mxArray *init)
{
    double seed = 0;

    if (mxIsNumeric(init) == 0
    ||  mxIsComplex(init) == 1
    ||  mxIsScalar(init)  == 0)
        mexErrMsgTxt("SEED must be a real, positive, integer scalar.");

    seed = mxGetScalar(init);

    if (seed < 1 || is_integer_64(seed) == 0)
        mexErrMsgTxt("SEED must be a real, positive, integer scalar.");

    if (seed > (double)0xFFFFFFFF)
        mexErrMsgTxt("SEED must be less than 2^32.");

    return seed;
}


void mexFunction(int            n_out,
                 mxArray       *out[],
                 int            n_in,
                 const mxArray *in[])
{
    double   *seq  = NULL;
    size_t    len  = 0;
    uint32_t  seed = 0;

    narginchk(n_in, 2, 2);
    nargoutchk(n_out, 0, 1);

    out[0] = allocate_output(in[0]);
    seq    = mxGetPr(out[0]);
    len    = mxGetNumberOfElements(out[0]);
    seed   = read_seed(in[1]);

    ordpat_xorshift_rand_double(seq, len, seed);
}
