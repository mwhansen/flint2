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

int
main(void)
{
    int i, j;
    FLINT_TEST_INIT(state);

    flint_printf("minpoly....");
    fflush(stdout);

    /* Test minpoly over fields */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fq_ctx_t ctx;
        fq_t sigma, e, t;
        fmpz_mod_poly_t mpoly;

        fq_ctx_randtest(ctx, state);

        fmpz_mod_poly_init(mpoly, fq_ctx_prime(ctx));

        fq_init(sigma, ctx);
        fq_init(e, ctx);
        fq_init(t, ctx);
        fq_randtest_not_zero(sigma, state, ctx);

        fq_minpoly(mpoly, sigma, ctx);

        fq_set_fmpz(e, mpoly->coeffs + mpoly->length - 1, ctx);
        for (j = mpoly->length - 2; j >= 0; j--)
        {
            fq_mul(e, e, sigma, ctx);
            fq_set_fmpz(t, mpoly->coeffs + j, ctx);
            fq_add(e, e, t, ctx);
        }
        
        if (!fq_is_zero(e, ctx))
        {
            flint_printf("FAIL:\n");
            fq_ctx_print(ctx); flint_printf("\n\n");
            fq_print_pretty(sigma, ctx); flint_printf("\n\n");
            fmpz_mod_poly_print_pretty(mpoly, "x"), flint_printf("\n\n");
            fq_print_pretty(e, ctx), flint_printf("\n\n");
            abort();
        }

        fmpz_mod_poly_clear(mpoly);
        fq_clear(sigma, ctx);
        fq_clear(e, ctx);
        fq_clear(t, ctx);
        fq_ctx_clear(ctx);
    }

    /* Test minpoly over non-fields */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fq_ctx_t ctx;
        fmpz_t p;
        fq_t sigma, e, t;
        fmpz_mod_poly_t modulus, mpoly;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        fmpz_mod_poly_init(modulus, p);
        do
        {
            fmpz_mod_poly_randtest(modulus, state, n_randint(state, 30));
        } while (modulus->length < 2);
        
        fq_ctx_init_modulus(ctx, modulus, "a");

        fmpz_mod_poly_init(mpoly, fq_ctx_prime(ctx));

        fq_init(sigma, ctx);
        fq_init(e, ctx);
        fq_init(t, ctx);
        fq_randtest_not_zero(sigma, state, ctx);

        fq_minpoly(mpoly, sigma, ctx);

        fq_set_fmpz(e, mpoly->coeffs + mpoly->length - 1, ctx);
        for (j = mpoly->length - 2; j >= 0; j--)
        {
            fq_mul(e, e, sigma, ctx);
            fq_set_fmpz(t, mpoly->coeffs + j, ctx);
            fq_add(e, e, t, ctx);
        }
        
        if (!fq_is_zero(e, ctx))
        {
            flint_printf("FAIL:\n");
            fq_ctx_print(ctx); flint_printf("\n\n");
            fq_print_pretty(sigma, ctx); flint_printf("\n\n");
            fmpz_mod_poly_print_pretty(mpoly, "x"), flint_printf("\n\n");
            fq_print_pretty(e, ctx), flint_printf("\n\n");
            abort();
        }

        fmpz_mod_poly_clear(mpoly);
        fmpz_mod_poly_clear(modulus);
        fmpz_clear(p);
        fq_clear(sigma, ctx);
        fq_clear(e, ctx);
        fq_clear(t, ctx);
        fq_ctx_clear(ctx);
    }
    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

