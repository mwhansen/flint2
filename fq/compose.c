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

#include "fq.h"

void
fq_compose(fq_t rop, const fq_t op, const fq_t xp, const fq_ctx_t ctx)
{
    fmpz* xpp;
    slong d = fq_ctx_degree(ctx);
    
    if (fq_is_zero(op, ctx))
    {
        fq_zero(rop, ctx);
        return;
    }

    if (op->length == 1)
    {
        fq_set(rop, op, ctx);
        return;
    }

    if (rop == op)
    {
        fq_t tmp;
        fq_init(tmp, ctx);
        fq_compose(tmp, op, xp, ctx);
        fq_swap(rop, tmp, ctx);
        fq_clear(tmp, ctx);
        return;
    }

    xpp = _fmpz_vec_init(d);
    _fmpz_vec_set(xpp, xp->coeffs, xp->length);
    _fmpz_vec_zero(xp->coeffs + xp->length, d - xp->length);
    
    fmpz_poly_fit_length(rop, d);
    _fmpz_mod_poly_compose_mod_brent_kung_preinv(rop->coeffs,
             op->coeffs, op->length, xpp,
             ctx->modulus->coeffs, ctx->modulus->length,
             ctx->inv->coeffs, ctx->inv->length,
             fq_ctx_prime(ctx));
    _fmpz_poly_set_length(rop, d);
    _fmpz_poly_normalise(rop);

    _fmpz_vec_clear(xpp, d);
}
