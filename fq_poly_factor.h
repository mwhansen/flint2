/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2013 Mike Hansen

******************************************************************************/

#ifndef FQ_POLY_FACTOR_H
#define FQ_POLY_FACTOR_H

static __inline__ int FQ_POLY_ITERATED_FROBENIUS_CUTOFF(const fq_ctx_t ctx, slong length)
{
    int result;
    fmpz_t q;
    fmpz_init(q);
    fq_ctx_order(q, ctx);
    if ( fmpz_sizeinbase(q, 2) < 3 * (n_sqrt(length) + 1))
        result = 1;
    else
        result = 0;
    fmpz_clear(q);
    return result;
}

/* Cutoff between von zur Gathen - Shoup and Kaltofen - Shoup */
static __inline__ int
FQ_POLY_FACTOR_EQUAL_DEG_USE_VZGS(const fq_poly_t poly,
                                  const fq_ctx_t ctx)
{
    double cutoff;
    cutoff = (0.088 * fmpz_bits(fq_ctx_prime(ctx)) +
              1.170 * FLINT_BIT_COUNT(fq_ctx_degree(ctx)) -
              1.840);
    if (sqrt((double) poly->length) > cutoff)
        return 1;
    else
        return 0;
}


#ifdef T
#undef T
#endif

#define T fq
#define CAP_T FQ
#include "fq_poly_factor_templates.h"
#undef CAP_T
#undef T

#endif
