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

#include "nmod_poly.h"

int
main(void)
{
    int i, j, k, result;
    FLINT_TEST_INIT(state);

    flint_printf("berlekamp_massey....");
    fflush(stdout);

    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        mp_limb_t p, lcinv;
        nmod_poly_t poly, res, q, r;
        mp_limb_t *c;
        slong d;

        p = n_randtest_prime(state, 0);

        nmod_poly_init(res, p);
        nmod_poly_init(poly, p);
        nmod_poly_init(q, p);
        nmod_poly_init(r, p);

        /* Generate random polynomial */
        nmod_poly_randtest_monic(poly, state, n_randint(state, 50) + 2);

        lcinv = n_invmod(poly->coeffs[poly->length - 1], p);

        d = poly->length - 1;

        /* Generate linearly recurrent sequence */
        c = _nmod_vec_init(2 * d);
        for (j = 0; j < d; j++)
        {
            c[j] = n_randint(state, p - 1) + 1;
        }
        for (j = 0; j < d; j++)
        {
            c[d + j] = 0;
            for (k = 0; k < d; k++)
            {
                c[d + j] = n_addmod(c[d + j],
                                    n_mulmod2_preinv(poly->coeffs[k], c[j + k],
                                                     poly->mod.n,
                                                     poly->mod.ninv),
                                    p);
            }
            c[d + j] = n_mulmod2_preinv(c[d + j], lcinv, poly->mod.n,
                                        poly->mod.ninv);
            c[d + j] = n_negmod(c[d + j], p);
        }

        nmod_poly_berlekamp_massey(res, c, d);
        nmod_poly_divrem(q, r, poly, res);

        result = (nmod_poly_is_zero(r));
        if (nmod_poly_is_irreducible(poly))
        {
            nmod_poly_make_monic(res, res);
            result = result && nmod_poly_equal(res, poly);
        }

        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("p = %wd\n", p);
            nmod_poly_print(poly); flint_printf("\n\n");
            nmod_poly_print(res); flint_printf("\n\n");
            abort();
        }

        _nmod_vec_clear(c);
        nmod_poly_clear(res);
        nmod_poly_clear(poly);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

