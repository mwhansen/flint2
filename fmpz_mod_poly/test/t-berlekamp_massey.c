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

#include "fmpz_mod_poly.h"

int
main(void)
{
    int i, j, k, result;
    FLINT_TEST_INIT(state);

    flint_printf("berlekamp_massey....");
    fflush(stdout);

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t p, lcinv;
        fmpz_mod_poly_t poly, res;
        fmpz *c;
        slong d;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        fmpz_set_ui(p, 5);

        fmpz_mod_poly_init(res, p);
        fmpz_mod_poly_init(poly, p);

        /* Generate random polynomial */
        fmpz_mod_poly_randtest_irreducible(poly, state, n_randint(state, 20) + 2);

        fmpz_init(lcinv);
        fmpz_invmod(lcinv, poly->coeffs + poly->length - 1, p);

        d = poly->length - 1;

        /* Generate linearly recurrent sequence */
        c = _fmpz_vec_init(2 * d);
        for (j = 0; j < d; j++)
        {
            fmpz_randm(c + j, state, p);
            if (fmpz_is_zero(c + j))
                fmpz_one(c + j);
        }
        for (j = 0; j < d; j++)
        {
            fmpz_zero(c + d + j);
            for (k = 0; k < d; k++)
            {
                fmpz_addmul(c + j + d,
                            poly->coeffs + k,
                            c + j + k);
            }
            fmpz_mul(c + j + d, c + j + d, lcinv);
            fmpz_neg(c + j + d, c + j + d);
            fmpz_mod(c + j + d, c + j + d, p);
        }

        fmpz_mod_poly_berlekamp_massey(res, c, d);
        fmpz_mod_poly_make_monic(res, res);
        fmpz_mod_poly_scalar_mul_fmpz(res, res, poly->coeffs + poly->length - 1);

        result = (fmpz_mod_poly_equal(res, poly));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("p = "); fmpz_print(p); flint_printf("\n\n");
            fmpz_mod_poly_print_pretty(poly, "x"), flint_printf("\n\n");
            _fmpz_vec_print(c, 2 * d); flint_printf("\n\n");
            fmpz_mod_poly_print_pretty(res, "x"), flint_printf("\n\n");
            abort();
        }

        _fmpz_vec_clear(c, 2 * d);
        fmpz_mod_poly_clear(res);
        fmpz_mod_poly_clear(poly);
        fmpz_clear(p);
        fmpz_clear(lcinv);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

