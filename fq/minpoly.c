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
fq_minpoly(fmpz_mod_poly_t g, const fq_t op, const fq_ctx_t ctx)
{
    int i, j, nonzero;
    flint_rand_t state;
    fmpz *c, *v, *op_coeffs;
    fmpz_mod_poly_t tau, g_prime, temp;
    ulong n = fq_ctx_degree(ctx), k;
    
    flint_randinit(state);
    
    fmpz_mod_poly_zero(g);
    fmpz_mod_poly_set_coeff_ui(g, 0, 1);
    
    fmpz_mod_poly_init(tau, fq_ctx_prime(ctx));
    fmpz_mod_poly_set_coeff_ui(tau, 0, 1);

    fmpz_mod_poly_init(temp, fq_ctx_prime(ctx));
    fmpz_mod_poly_init(g_prime, fq_ctx_prime(ctx));

    c = _fmpz_vec_init(2 * n);
    v = _fmpz_vec_init(n);
    op_coeffs = _fmpz_vec_init(n);
    _fmpz_vec_set(op_coeffs, op->coeffs, op->length);

    /* Use (1, 0, ...) for the first iteration */
    fmpz_one(v);

    while (!fmpz_mod_poly_is_zero(tau))
    {
        k = 2*(n - fmpz_mod_poly_degree(g));
        fq_power_projection(c, op, v, k, ctx);

        fmpz_mod_poly_berlekamp_massey(g_prime, c, k / 2);

        /* g = g * g_prime */
        fmpz_mod_poly_mul(temp, g, g_prime);
        fmpz_mod_poly_swap(g, temp);

        if (fmpz_mod_poly_degree(g) == n)
            break;

        /* tau = tau * g_prime(op) */
        if (g_prime->length > 1)
        {
            fmpz_mod_poly_fit_length(temp, n);
            for (j = op->length; j < n; j++)
                fmpz_zero(op->coeffs + j);
            _fmpz_mod_poly_compose_mod_brent_kung_preinv(
                temp->coeffs,
                g_prime->coeffs, g_prime->length,
                op_coeffs,
                ctx->modulus->coeffs, ctx->modulus->length,
                ctx->inv->coeffs, ctx->inv->length,
                fq_ctx_prime(ctx));
            _fmpz_mod_poly_set_length(temp, n);
            _fmpz_mod_poly_normalise(temp);
            fmpz_mod_poly_swap(g_prime, temp);
        }
        fmpz_mod_poly_mulmod_preinv(temp, tau, g_prime, ctx->modulus, ctx->inv);
        fmpz_mod_poly_swap(tau, temp);

        /* Choose v \in K^n at random for next iteration */
        nonzero = 0;
        while (nonzero == 0)
        {
            nonzero = 0;
            for (i = 0; i < n; i++)
            {
                fmpz_randm(v + i, state, fq_ctx_prime(ctx));
                nonzero = nonzero || !fmpz_is_zero(v + i);
            }
        }
        fmpz_mod_poly_mul_transposed_preinv(v, tau, v, ctx->modulus, ctx->inv);
    }

    _fmpz_vec_clear(c, 2 * n);
    _fmpz_vec_clear(v, n);
    _fmpz_vec_clear(op_coeffs, n);

    fmpz_mod_poly_clear(temp);
    fmpz_mod_poly_clear(g_prime);
    fmpz_mod_poly_clear(tau);
    flint_randclear(state);
    

}
