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

    Copyright (C) 2014 Mike Hansen

******************************************************************************/


#ifdef T

#include "templates.h"

ulong
TEMPLATE(T, ctx_randtest_with_e)(TEMPLATE(T, ctx_t) ctx,
                                 flint_rand_t state)
{
    ulong e;
    fmpz_t qm1;
    slong i, count = 0;
    fmpz_factor_t factor;

    TEMPLATE(T, ctx_randtest)(ctx, state);
    while (TEMPLATE(T, ctx_degree)(ctx) == 1)
    {
        TEMPLATE(T, ctx_clear)(ctx);
        TEMPLATE(T, ctx_randtest)(ctx, state);
    }

    /* Choose an e that divides q - 1*/
    fmpz_init(qm1);
    fmpz_factor_init(factor);
    TEMPLATE(T, ctx_order)(qm1, ctx);
    fmpz_sub_ui(qm1, qm1, 1);
    fmpz_factor(factor, qm1);

    /* Pick a random prime factor less than 32 */
    for (i = 0; i < factor->num; i++)
    {
        if (fmpz_cmp_ui(factor->p + i, 32) <= 0)
        {
            count += 1;
        }
    }
    if (count == 0)
    {
        e = fmpz_get_ui(factor->p);
    }
    else
    {
        e = fmpz_get_ui(factor->p + n_randint(state, count));
    }
    fmpz_factor_clear(factor);
    fmpz_clear(qm1);
    return e;
}

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("root_ds... ");
    fflush(stdout);

    /* Test root_ds_s_ell */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        ulong t;
        fmpz_t s, ell, e, n, tt;
        TEMPLATE(T, ctx_t) ctx;
        
        t = TEMPLATE(T, ctx_randtest_with_e)(ctx, state);

        fmpz_init(s);
        fmpz_init(ell);
        fmpz_init(e);
        fmpz_init(n);
        fmpz_init_set_ui(tt, t);

        __TEMPLATE(T, root_ds_s_ell)(s, ell, e, t, ctx);

        result = 1;

        fmpz_mul(n, s, ell);
        result = result && (fmpz_cmp_ui(n, TEMPLATE(T, ctx_degree)(ctx)) == 0);

        fmpz_mod(n, TEMPLATE(T, ctx_prime)(ctx), tt);
        fmpz_powm(n, n, s, tt);
        result = result && fmpz_is_one(n);

        /* Test e = (p^s - 1) / t */
        fmpz_pow_ui(n, TEMPLATE(T, ctx_prime)(ctx), fmpz_get_ui(s));
        fmpz_sub_ui(n, n, 1);
        fmpz_divexact_ui(n, n, t);
        result = result && fmpz_equal(n, e);

        if (!result)
        {
            flint_printf("FAIL (s_ell):\n\n");
            TEMPLATE(T, ctx_print)(ctx);
            flint_printf("\nt = %wu\n", t);
            flint_printf("s = "), fmpz_print(s), flint_printf("\n");
            flint_printf("ell = "), fmpz_print(ell), flint_printf("\n");
            flint_printf("e = "), fmpz_print(e), flint_printf("\n");
            abort();
        }

        fmpz_clear(s);
        fmpz_clear(ell);
        fmpz_clear(e);
        fmpz_clear(n);

        TEMPLATE(T, ctx_clear)(ctx);
    }


    /* Test root_ds_xi1zeta1 */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        int j;
        fmpz_t s;
        TEMPLATE(T, t) a, b, xi1, zeta1, xp, x;
        TEMPLATE(T, ctx_t) ctx;
        
        TEMPLATE(T, ctx_randtest)(ctx, state);

        fmpz_init(s);

        fmpz_randtest_unsigned(s, state, 5);
        fmpz_add_ui(s, s, 1);
        TEMPLATE(T, init)(a, ctx);
        TEMPLATE(T, init)(b, ctx);
        TEMPLATE(T, init)(x, ctx);
        TEMPLATE(T, init)(xi1, ctx);
        TEMPLATE(T, init)(zeta1, ctx);
        TEMPLATE(T, init)(xp, ctx);
        TEMPLATE(T, randtest)(a, state, ctx);

        __TEMPLATE(T, root_ds_xi1zeta1)(xi1, zeta1, xp, a, s, ctx);

        result = 1;

        /* Test xp == x ^ p */
        TEMPLATE(T, gen)(x, ctx);
        TEMPLATE(T, pow)(b, x, TEMPLATE(T, ctx_prime)(ctx), ctx);
        result = result && TEMPLATE(T, equal)(b, xp, ctx);

        /* Test xi1 = x ^ (p^s) */
        for (j = 1; j < fmpz_get_ui(s); j++)
        {
            TEMPLATE(T, pow)(b, b, TEMPLATE(T, ctx_prime)(ctx), ctx);
        }
        result = result && TEMPLATE(T, equal)(b, xi1, ctx);

        /* Test zeta1 = a ^ (p^s) */
        TEMPLATE(T, set)(b, a, ctx);
        for (j = 0; j < fmpz_get_ui(s); j++)
        {
            TEMPLATE(T, pow)(b, b, TEMPLATE(T, ctx_prime)(ctx), ctx);
        }
        result = result && TEMPLATE(T, equal)(b, zeta1, ctx);

        if (!result)
        {
            flint_printf("FAIL (xi1zeta1):\n\n");
            TEMPLATE(T, ctx_print)(ctx);
            flint_printf("\ns = %wu\n", fmpz_get_ui(s));;
            flint_printf("a = "), TEMPLATE(T, print_pretty)(a, ctx), flint_printf("\n");
            flint_printf("xi1 = "), TEMPLATE(T, print_pretty)(xi1, ctx), flint_printf("\n");
            flint_printf("zeta1 = "), TEMPLATE(T, print_pretty)(zeta1, ctx), flint_printf("\n");
            flint_printf("xp = "), TEMPLATE(T, print_pretty)(xp, ctx), flint_printf("\n");
            abort();
        }

        TEMPLATE(T, clear)(a, ctx);
        TEMPLATE(T, clear)(b, ctx);
        TEMPLATE(T, clear)(x, ctx);
        TEMPLATE(T, clear)(xp, ctx);
        TEMPLATE(T, clear)(xi1, ctx);
        TEMPLATE(T, clear)(zeta1, ctx);
        
        fmpz_clear(s);

        TEMPLATE(T, ctx_clear)(ctx);
    }
    
    /* Test root_ds_xizeta */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        slong ii;
        int j;
        fmpz_t s, n, temp;
        TEMPLATE(T, t) a, b, xi, zeta, xi1, zeta1, xp, x;
        TEMPLATE(T, ctx_t) ctx;
        
        TEMPLATE(T, ctx_randtest)(ctx, state);

        fmpz_init(s);

        fmpz_randtest_unsigned(s, state, 5);
        fmpz_add_ui(s, s, 1);
        TEMPLATE(T, init)(a, ctx);
        TEMPLATE(T, init)(b, ctx);
        TEMPLATE(T, init)(x, ctx);
        TEMPLATE(T, init)(xi, ctx);
        TEMPLATE(T, init)(zeta, ctx);
        TEMPLATE(T, init)(xi1, ctx);
        TEMPLATE(T, init)(zeta1, ctx);
        TEMPLATE(T, init)(xp, ctx);
        TEMPLATE(T, randtest)(a, state, ctx);

        ii = n_randint(state, 10) + 1;
        __TEMPLATE(T, root_ds_xi1zeta1)(xi1, zeta1, xp, a, s, ctx);
        __TEMPLATE(T, root_ds_xizeta)(xi, zeta, ii, xi1, zeta1, ctx);

        result = 1;

        /* Test xi == x ^ (p^(ii*s) */
        fmpz_init(n);
        TEMPLATE(T, gen)(x, ctx);
        fmpz_pow_ui(n, TEMPLATE(T, ctx_prime)(ctx), ii * fmpz_get_ui(s));
        TEMPLATE(T, pow)(b, x, n, ctx);
        result = result && TEMPLATE(T, equal)(b, xi, ctx);

        /* Testa zeta = a ^ (p^s + p^(2*s) + ... + p^(i*s) */
        fmpz_init(temp);
        fmpz_zero(n);
        for (j = 1; j <= ii; j++)
        {
            fmpz_pow_ui(temp, TEMPLATE(T, ctx_prime)(ctx), j * fmpz_get_ui(s));
            fmpz_add(n, n, temp);
        }
        TEMPLATE(T, pow)(b, a, n, ctx);
        result = result && TEMPLATE(T, equal)(b, zeta, ctx);

        if (!result)
        {
            flint_printf("FAIL (xizeta):\n\n");
            TEMPLATE(T, ctx_print)(ctx);
            flint_printf("\ns = %wu\n", fmpz_get_ui(s));
            flint_printf("\ni = %wd\n", ii);
            flint_printf("a = "), TEMPLATE(T, print_pretty)(a, ctx), flint_printf("\n");
            flint_printf("xi = "), TEMPLATE(T, print_pretty)(xi, ctx), flint_printf("\n");
            flint_printf("zeta = "), TEMPLATE(T, print_pretty)(zeta, ctx), flint_printf("\n");
            abort();
        }

        TEMPLATE(T, clear)(a, ctx);
        TEMPLATE(T, clear)(b, ctx);
        TEMPLATE(T, clear)(x, ctx);
        TEMPLATE(T, clear)(xp, ctx);
        TEMPLATE(T, clear)(xi1, ctx);
        TEMPLATE(T, clear)(zeta1, ctx);
        TEMPLATE(T, clear)(xi, ctx);
        TEMPLATE(T, clear)(zeta, ctx);
        
        fmpz_clear(s);
        fmpz_clear(n);
        fmpz_clear(temp);

        TEMPLATE(T, ctx_clear)(ctx);
    }


    /* Test root_ds_qm1ot */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        ulong t;
        fmpz_t s, ell, e, n, tt;
        TEMPLATE(T, t) a, b, c;
        TEMPLATE(T, ctx_t) ctx;
        
        t = TEMPLATE(T, ctx_randtest_with_e)(ctx, state);

        fmpz_init(s);
        fmpz_init(ell);
        fmpz_init(e);
        fmpz_init(n);
        fmpz_init_set_ui(tt, t);

        TEMPLATE(T, init)(a, ctx);
        TEMPLATE(T, init)(b, ctx);
        TEMPLATE(T, init)(c, ctx);
        TEMPLATE(T, randtest)(a, state, ctx);

        __TEMPLATE(T, root_ds_s_ell)(s, ell, e, t, ctx);
        __TEMPLATE(T, root_ds_qm1ot)(b, a, s, ell, e, t, ctx);

        TEMPLATE(T, ctx_order)(n, ctx);
        fmpz_sub_ui(n, n, 1);
        fmpz_divexact_ui(n, n, t);

        TEMPLATE(T, pow)(c, a, n, ctx);
        
        result = TEMPLATE(T, equal)(b, c, ctx);

        if (!result)
        {
            flint_printf("FAIL (qm1ot):\n\n");
            TEMPLATE(T, ctx_print)(ctx);
            flint_printf("\nt = %wu\n", t);
            flint_printf("n = "), fmpz_print(n), flint_printf("\n");
            flint_printf("a = "), TEMPLATE(T, print_pretty)(a, ctx), flint_printf("\n");
            flint_printf("b = "), TEMPLATE(T, print_pretty)(b, ctx), flint_printf("\n");
            flint_printf("c = "), TEMPLATE(T, print_pretty)(c, ctx), flint_printf("\n");
            abort();
        }

        TEMPLATE(T, clear)(a, ctx);
        TEMPLATE(T, clear)(b, ctx);
        TEMPLATE(T, clear)(c, ctx);
        
        fmpz_clear(s);
        fmpz_clear(ell);
        fmpz_clear(e);
        fmpz_clear(n);

        TEMPLATE(T, ctx_clear)(ctx);
    }
    

    /* Test on t-th powers */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        ulong e;
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, t) a, b, c;

        e = TEMPLATE(T, ctx_randtest_with_e)(ctx, state);
        
        TEMPLATE(T, init)(a, ctx);
        TEMPLATE(T, init)(b, ctx);
        TEMPLATE(T, init)(c, ctx);

        TEMPLATE(T, randtest_not_zero)(a, state, ctx);

        TEMPLATE(T, pow_ui)(a, a, e, ctx);

        result = TEMPLATE(T, root_ds)(b, a, e, ctx);

        TEMPLATE(T, pow_ui)(c, b, e, ctx);

        result = result && TEMPLATE(T, equal)(a, c, ctx);
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            TEMPLATE(T, ctx_print)(ctx);
            flint_printf("\ne = %wu\n", e);
            flint_printf("a = "), TEMPLATE(T, print_pretty)(a, ctx), flint_printf("\n");
            flint_printf("b = "), TEMPLATE(T, print_pretty)(b, ctx), flint_printf("\n");
            flint_printf("c = "), TEMPLATE(T, print_pretty)(c, ctx), flint_printf("\n");
            abort();
        }

        TEMPLATE(T, clear)(a, ctx);
        TEMPLATE(T, clear)(b, ctx);
        TEMPLATE(T, clear)(c, ctx);

        TEMPLATE(T, ctx_clear)(ctx);
    }

    /* Test on random elements */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        ulong e;
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, t) a, b, c;

        e = TEMPLATE(T, ctx_randtest_with_e)(ctx, state);
        
        TEMPLATE(T, init)(a, ctx);
        TEMPLATE(T, init)(b, ctx);
        TEMPLATE(T, init)(c, ctx);

        TEMPLATE(T, randtest_not_zero)(a, state, ctx);

        result = TEMPLATE(T, root_ds)(b, a, e, ctx);
        
        if (result)
        {
            TEMPLATE(T, pow_ui)(c, b, e, ctx);
            result = TEMPLATE(T, equal)(a, c, ctx);
        }
        else
        {
            result = 1;
        }
        
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), TEMPLATE(T, print_pretty)(a, ctx), flint_printf("\n");
            flint_printf("b = "), TEMPLATE(T, print_pretty)(b, ctx), flint_printf("\n");
            flint_printf("c = "), TEMPLATE(T, print_pretty)(c, ctx), flint_printf("\n");
            abort();
        }

        TEMPLATE(T, clear)(a, ctx);
        TEMPLATE(T, clear)(b, ctx);
        TEMPLATE(T, clear)(c, ctx);

        TEMPLATE(T, ctx_clear)(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}


#endif
