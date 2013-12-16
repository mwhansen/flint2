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

#include "fq_nmod.h"
#include "fq_nmod_vec.h"

void
fq_nmod_power_projection(mp_ptr c, const fq_nmod_t sigma, mp_srcptr v, slong l,
                         const fq_nmod_ctx_t ctx)
{
    int nlimbs;
    mp_limb_t *vv;
    fq_nmod_struct *sigma_j;
    slong i, j, counter, k = n_sqrt(l), k_prime = (1 + ((l - 1) / k));
    slong n = fq_nmod_ctx_degree(ctx);

    /* Handle aliasing */
    vv = _nmod_vec_init(n);
    _nmod_vec_zero(vv, n);
    _nmod_vec_set(vv, v, n);
    
    /* Compute 1, sigma, sigma^2, ..., sigma^k */
    sigma_j = _fq_nmod_vec_init2(k + 1, ctx);
    fq_nmod_one(sigma_j + 0, ctx);
    fq_nmod_set(sigma_j + 1, sigma, ctx);
    for (i = 2; i <= k; i++)
    {
        fq_nmod_mul(sigma_j + i, sigma_j + (i - 1), sigma_j + 1, ctx);
    }

    /* Compute power projection coefficients */
    counter = 0;
    nlimbs = _nmod_vec_dot_bound_limbs(n, ctx->mod);
    for (i = 0; i < k_prime; i++)
    {
        for (j = 0; j < k; j++)
        {
            if (counter >= l)
                goto exit;
            c[counter] = _nmod_vec_dot(vv, (sigma_j + j)->coeffs,
                                       FLINT_MIN(n, (sigma_j + j)->length),
                                       ctx->mod, nlimbs);
            counter += 1;
        }
        _nmod_poly_mulmod_transposed_preinv(
            vv, (sigma_j + k)->coeffs, (sigma_j + k)->length,
            vv, ctx->modulus->coeffs, ctx->modulus->length,
            ctx->inv->coeffs, ctx->inv->length, ctx->mod);
    }

exit:
    _nmod_vec_clear(vv);
    _fq_nmod_vec_clear(sigma_j, k + 1, ctx);
}
