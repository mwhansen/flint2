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

    Copyright (C) 2012 Andres Goens
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

******************************************************************************/

#ifndef FQ_NMOD_POLY_H
#define FQ_NMOD_POLY_H

#include "fq_nmod.h"
#include "fq_nmod_mat.h"
#include "math.h"

#define FQ_NMOD_POLY_DIVREM_DIVCONQUER_CUTOFF  16
#define FQ_NMOD_COMPOSE_MOD_LENH_CUTOFF 6
#define FQ_NMOD_COMPOSE_MOD_PREINV_LENH_CUTOFF 6
#define FQ_NMOD_MUL_CLASSICAL_CUTOFF 6
#define FQ_NMOD_SQR_CLASSICAL_CUTOFF 6
#define FQ_NMOD_MULLOW_CLASSICAL_CUTOFF 6

#define FQ_NMOD_POLY_HGCD_CUTOFF 25
#define FQ_NMOD_POLY_SMALL_GCD_CUTOFF 110
#define FQ_NMOD_POLY_GCD_CUTOFF 120


#ifdef T
#undef T
#endif

#define T fq_nmod
#define CAP_T FQ_NMOD
#include "fq_poly_templates.h"
#undef CAP_T
#undef T

#include "fq_nmod_poly_factor.h"

static __inline__ int
FQ_NMOD_POLY_POWMOD_XQ_PREINV_USE_DIRECT(const fq_nmod_poly_t poly,
                                         const fq_nmod_ctx_t ctx)
{
    double cutoff;
    cutoff = (0.856 * fmpz_bits(fq_nmod_ctx_prime(ctx)) +
              3.795 * FLINT_BIT_COUNT(fq_nmod_ctx_degree(ctx)) -
              11.044);
    if (sqrt((double) poly->length) > cutoff)
        return 1;
    else
        return 0;
}

#endif
