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

void
fq_nmod_minpoly(nmod_poly_t g, const fq_nmod_t op, const fq_nmod_ctx_t ctx)
{
    int i, j, nonzero;
    flint_rand_t state;
    mp_limb_t *c, *v, *op_coeffs;
    nmod_poly_t tau, g_prime, temp;
    slong n = fq_nmod_ctx_degree(ctx), k;
    
    flint_randinit(state);
    
    nmod_poly_zero(g);
    nmod_poly_set_coeff_ui(g, 0, 1);
    
    nmod_poly_init(tau, ctx->mod.n);
    nmod_poly_set_coeff_ui(tau, 0, 1);

    nmod_poly_init(temp, ctx->mod.n);
    nmod_poly_init(g_prime, ctx->mod.n);

    c = _nmod_vec_init(2 * n);
    v = _nmod_vec_init(n);
    _nmod_vec_zero(v, n);
    op_coeffs = _nmod_vec_init(n);
    _nmod_vec_zero(op_coeffs, n);
    _nmod_vec_set(op_coeffs, op->coeffs, op->length);

    /* Use (1, 0, ...) for the first iteration */
    v[0] = 1;

    while (!nmod_poly_is_zero(tau))
    {
        k = 2*(n - nmod_poly_degree(g));
        fq_nmod_power_projection(c, op, v, k, ctx);

        nmod_poly_berlekamp_massey(g_prime, c, k / 2);

        /* g = g * g_prime */
        nmod_poly_mul(temp, g, g_prime);
        nmod_poly_swap(g, temp);

        if (nmod_poly_degree(g) == n)
            break;

        /* tau = tau * g_prime(op) */
        if (g_prime->length > 1)
        {
            nmod_poly_fit_length(temp, n);
            for (j = op->length; j < n; j++)
                op->coeffs[j] = 0;
            _nmod_poly_compose_mod_brent_kung_preinv(
                temp->coeffs,
                g_prime->coeffs, g_prime->length,
                op_coeffs,
                ctx->modulus->coeffs, ctx->modulus->length,
                ctx->inv->coeffs, ctx->inv->length,
                ctx->mod);
            _nmod_poly_set_length(temp, n);
            _nmod_poly_normalise(temp);
            nmod_poly_swap(g_prime, temp);
        }
        nmod_poly_mulmod_preinv(temp, tau, g_prime,
                                ctx->modulus, ctx->inv);
        nmod_poly_swap(tau, temp);

        /* Choose v \in K^n at random for next iteration */
        nonzero = 0;
        while (!nonzero)
        {
            nonzero = 0;
            for (i = 0; i < n; i++)
            {
                v[i] = n_randint(state, ctx->mod.n);
                nonzero = nonzero || v[i] != 0;
            }
        }
        nmod_poly_mulmod_transposed_preinv(v, tau, v,
                                           ctx->modulus, ctx->inv);
    }

    _nmod_vec_clear(c);
    _nmod_vec_clear(v);
    _nmod_vec_clear(op_coeffs);

    nmod_poly_clear(temp);
    nmod_poly_clear(g_prime);
    nmod_poly_clear(tau);
    flint_randclear(state);
    
}
