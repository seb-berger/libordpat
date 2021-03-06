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

#include <stdint.h>
#include <mex.h>

#include "ordpat.h"
#include "mex_helpers.h"


/* MEX function: TABLE = CREATE_LOOKUP_TABLE_C(ORD) */

void mexFunction(int            n_out,
                 mxArray       *out[],
                 int            n_in,
                 const mxArray *in[])
{
    uint64_t     *tab;
    size_t        len;
    unsigned int  ord;
    int           ret;

    narginchk(n_in, 1, 1);
    nargoutchk(n_out, 0, 1);

    ord = mx_to_pattern_order(in[0]);
    if (ord > 10)
        mexErrMsgTxt("ORD is limited to a maximum of 10.");

    len    = ordpat_factorial(ord);
    out[0] = mxCreateNumericMatrix(ord, len, mxUINT64_CLASS, mxREAL);
    tab    = mxGetData(out[0]);

    ret = ordpat_create_lookup_table(tab, len * ord, ord);

    if (ret)
        catch_libordpat_error(ret);
}
