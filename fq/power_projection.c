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
#include "fq_vec.h"

static __inline__ void
__fmpz_vec_dot(fmpz_t c, const fmpz *v, const fmpz *sigma, ulong n, const fmpz_t p)
{
    int i;
    fmpz_zero(c);
    for (i = 0; i < n; i++)
        fmpz_addmul(c, v + i, sigma + i);
    fmpz_mod(c, c, p);
}

void
fq_power_projection(fmpz *c, const fq_t sigma, const fmpz *v, ulong l,
                    const fq_ctx_t ctx)
{
    fmpz *vv;
    fq_struct *sigma_j;
    ulong i, j, k = n_sqrt(l), k_prime = (1 + ((l - 1) / k));
    ulong n = fq_ctx_degree(ctx);

    /* Handle aliasing */
    vv = _fmpz_vec_init(n);
    _fmpz_vec_set(vv, v, n);
    
    /* Compute 1, sigma, sigma^2, ..., sigma^k */
    sigma_j = _fq_vec_init2(k + 1, ctx);
    fq_one(sigma_j + 0, ctx);
    fq_set(sigma_j + 1, sigma, ctx);
    for (i = 2; i <= k; i++)
    {
        fq_mul(sigma_j + i, sigma_j + (i - 1), sigma_j + 1, ctx);
    }

    /* Compute power projection coefficients */
    for (i = 0; i < k_prime; i++)
    {
        for (j = 0; j < k; j++)
        {
            if (i == k_prime - 1 && i*k + j >= l)
                goto exit;
            __fmpz_vec_dot(c + (i*k + j),
                           vv, (sigma_j + j)->coeffs,
                           FLINT_MIN(n, (sigma_j + j)->length),
                           fq_ctx_prime(ctx));
        }
        _fmpz_mod_poly_mulmod_transposed_preinv(
            vv, (sigma_j + k)->coeffs, (sigma_j + k)->length,
            vv, ctx->modulus->coeffs, ctx->modulus->length,
            ctx->inv->coeffs, ctx->inv->length, fq_ctx_prime(ctx));
    }

exit:
    _fmpz_vec_clear(vv, n);
    _fq_vec_clear(sigma_j, k + 1, ctx);
}
