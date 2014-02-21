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

    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

******************************************************************************/

#include "fmpz_mod_poly.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("resultant_euclidean....");
    fflush(stdout);

    /* Check res(f, g) == (-1)^(deg f deg g) res(g, f) */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t f, g;
        fmpz_t x, y, n;

        fmpz_init(x);
        fmpz_init(y);
        fmpz_init(n);
        
        fmpz_init_set_ui(n, n_randtest_prime(state, 0));

        fmpz_mod_poly_init(f, n);
        fmpz_mod_poly_init(g, n);
        
        fmpz_mod_poly_randtest(f, state, n_randint(state, 200));
        fmpz_mod_poly_randtest(g, state, n_randint(state, 200));

        fmpz_mod_poly_resultant_euclidean(x, f, g);
        fmpz_mod_poly_resultant_euclidean(y, g, f);

        if ((fmpz_mod_poly_degree(f) * fmpz_mod_poly_degree(g)) % 2)
            if (fmpz_cmp_ui(y, 0) > 0)
                fmpz_sub(y, n, y);

        result = fmpz_equal(x, y);
        if (!result)
        {
            flint_printf("FAIL (res(f, g) == (-1)^(deg f deg g) res(g, f)):\n");
            fmpz_mod_poly_print_pretty(f, "x"), flint_printf("\n\n");
            fmpz_mod_poly_print_pretty(g, "x"), flint_printf("\n\n");
            flint_printf("x = ", x); fmpz_print(x); flint_printf("\n");
            flint_printf("y = ", y); fmpz_print(y); flint_printf("\n");
            flint_printf("n = ", n); fmpz_print(n); flint_printf("\n");
            abort();
        }
        
        fmpz_mod_poly_clear(f);
        fmpz_mod_poly_clear(g);

        fmpz_clear(x);
        fmpz_clear(y);
        fmpz_clear(n);
    }

    /* Check res(f h, g) == res(f, g) res(h, g) */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t f, g, h;
        fmpz_t x, y, z, n;

        fmpz_init(x);
        fmpz_init(y);
        fmpz_init(z);
        fmpz_init(n);
        
        fmpz_init_set_ui(n, n_randtest_prime(state, 0));

        fmpz_mod_poly_init(f, n);
        fmpz_mod_poly_init(g, n);
        fmpz_mod_poly_init(h, n);
        
        fmpz_mod_poly_randtest(f, state, n_randint(state, 200));
        fmpz_mod_poly_randtest(g, state, n_randint(state, 200));
        fmpz_mod_poly_randtest(h, state, n_randint(state, 200));

        fmpz_mod_poly_resultant_euclidean(y, f, g);
        fmpz_mod_poly_resultant_euclidean(z, h, g);
        fmpz_mul(y, y, z);
        fmpz_mod(y, y, n);
        fmpz_mod_poly_mul(f, f, h);
        fmpz_mod_poly_resultant_euclidean(x, f, g);

        result = fmpz_equal(x, y);
        if (!result)
        {
            flint_printf("FAIL (res(f h, g) == res(f, g) res(h, g)):\n");
            fmpz_mod_poly_print_pretty(f, "x"), flint_printf("\n\n");
            fmpz_mod_poly_print_pretty(g, "x"), flint_printf("\n\n");
            flint_printf("x = "); fmpz_print(x); flint_printf("\n");
            flint_printf("y = "); fmpz_print(y); flint_printf("\n");
            flint_printf("z = "); fmpz_print(z); flint_printf("\n");
            flint_printf("n = "); fmpz_print(n); flint_printf("\n");
            abort();
        }
        
        fmpz_mod_poly_clear(f);
        fmpz_mod_poly_clear(g);
        fmpz_mod_poly_clear(h);

        fmpz_clear(x);
        fmpz_clear(y);
        fmpz_clear(z);
        fmpz_clear(n);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
