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

    flint_printf("root_edf... ");
    fflush(stdout);

    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        ulong e;
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, t) a, b, c;

        TEMPLATE(T, ctx_randtest)(ctx, state);
        
        TEMPLATE(T, init)(a, ctx);
        TEMPLATE(T, init)(b, ctx);
        TEMPLATE(T, init)(c, ctx);

        TEMPLATE(T, randtest_not_zero)(a, state, ctx);
        e = n_randint(state, 20) + 1;
        TEMPLATE(T, pow_ui)(a, a, e, ctx);

        result = TEMPLATE(T, root_edf)(b, a, e, ctx);
        TEMPLATE(T, pow_ui)(c, b, e, ctx);
        
        result = result && TEMPLATE(T, equal)(a, c, ctx);
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

    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        ulong e;
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, t) a, b, c;

        TEMPLATE(T, ctx_randtest)(ctx, state);
        
        TEMPLATE(T, init)(a, ctx);
        TEMPLATE(T, init)(b, ctx);
        TEMPLATE(T, init)(c, ctx);

        TEMPLATE(T, randtest_not_zero)(a, state, ctx);

        e = n_randint(state, 20) + 1;

        result = TEMPLATE(T, root_edf)(b, a, e, ctx);
        
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
