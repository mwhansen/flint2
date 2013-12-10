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


#ifdef T

#include "templates.h"


int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("trace_frob_preinv....");
    fflush(stdout);

    for (i = 0; i < 5 * flint_test_multiplier(); i++)
    {
        int j;
        ulong d;
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, poly_t) a, aq, res1, res2, f, finv;
        fmpz_t q;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        fmpz_init(q);
        TEMPLATE(T, ctx_order)(q, ctx);

        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (aq, ctx);
        TEMPLATE(T, poly_init) (f, ctx);
        TEMPLATE(T, poly_init) (res1, ctx);
        TEMPLATE(T, poly_init) (res2, ctx);
        TEMPLATE(T, poly_init) (finv, ctx);

        d = n_randint(state, 10);

        do
        {
            TEMPLATE(T, poly_randtest_monic) (f, state,
                                              n_randint(state, 20) + 2, ctx);
        } while (!TEMPLATE(T, poly_is_squarefree)(f, ctx));
        
        TEMPLATE(T, poly_reverse)(finv, f, f->length, ctx);
        TEMPLATE(T, poly_inv_series_newton)(finv, finv, f->length, ctx);
        TEMPLATE(T, poly_randtest)(a, state, f->length - 1, ctx);

        TEMPLATE(T, poly_powmod_xq_preinv)(aq, f, finv, ctx);
        TEMPLATE(T, poly_trace_frob_preinv) (res1, a, d, aq, f, finv, ctx);

        TEMPLATE(T, poly_set)(res2, a, ctx);
        TEMPLATE(T, poly_set)(aq, a, ctx);
        for (j = 1; j < d + 1; j++)
        {
            TEMPLATE(T, poly_powmod_fmpz_sliding_preinv)(aq, aq, q, 0, f, finv, ctx);
            TEMPLATE(T, poly_add)(res2, res2, aq, ctx);
        }

        result = (TEMPLATE(T, poly_equal) (res1, res2, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("d: %wu\n", d);
            flint_printf("a:\n");
            TEMPLATE(T, poly_print) (a, ctx), flint_printf("\n\n");
            flint_printf("f:\n");
            TEMPLATE(T, poly_print) (f, ctx), flint_printf("\n\n");
            flint_printf("res1:\n");
            TEMPLATE(T, poly_print) (res1, ctx), flint_printf("\n\n");
            flint_printf("res2:\n");
            TEMPLATE(T, poly_print) (res2, ctx), flint_printf("\n\n");
            abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (aq, ctx);
        TEMPLATE(T, poly_clear) (f, ctx);
        TEMPLATE(T, poly_clear) (res1, ctx);
        TEMPLATE(T, poly_clear) (res2, ctx);
        TEMPLATE(T, poly_clear) (finv, ctx);
        fmpz_clear(q);

        TEMPLATE(T, ctx_clear) (ctx);
    }


    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}

#endif
