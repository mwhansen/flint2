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
    Copyright (C) 2014 Mike Hansen

******************************************************************************/

#include "fmpz.h"


int main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("rootmod....");
    fflush(stdout);

    for (i = 0; i < 100 * flint_test_multiplier(); i++) /* Test random integers */
    {
        int ans;
        fmpz_t a, b, c, p;
        slong n;
        mp_limb_t prime;

        prime = n_randint(state, UWORD(1) << (FLINT_BITS - 1));
        prime = n_nextprime(prime, 1);

        n = n_randint(state, 32) + 2;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(p);

        fmpz_set_ui(p, prime);
        fmpz_randm(a, state, p);

        ans = fmpz_rootmod(b, a, n, p);

        fmpz_powm_ui(c, b, n, p);

        result = (ans == 0 || fmpz_equal(a, c));
        if (!result)
        {
            flint_printf("FAIL (random):\n");
            flint_printf("p = "), fmpz_print(p), flint_printf("\n");
            flint_printf("a = "), fmpz_print(a), flint_printf("\n");
            flint_printf("b = "), fmpz_print(b), flint_printf("\n");
            flint_printf("c = "), fmpz_print(c), flint_printf("\n");
            flint_printf("n = %wd\n", n);
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(p);
    }

    for (i = 0; i < 100 * flint_test_multiplier(); i++) /* Test random n^th powers */
    {
        int ans;
        fmpz_t a, b, c, d, p;
        mp_limb_t prime;
        slong n;

        prime = n_randint(state, UWORD(1) << (FLINT_BITS - 1));
        prime = n_nextprime(prime, 1);

        n = n_randint(state, 32) + 2;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(d);
        fmpz_init(p);

        fmpz_set_ui(p, prime);
        do 
            fmpz_randm(b, state, p);
        while (fmpz_is_zero(b));

        /* a = b ^ n mod p */
        fmpz_powm_ui(a, b, n, p);

        ans = fmpz_rootmod(c, a, n, p);

        fmpz_powm_ui(d, c, n, p);

        result = (ans && fmpz_equal(a, d));
        if (!result)
        {
            flint_printf("FAIL (squares):\n");
            flint_printf("p            = "), fmpz_print(p), flint_printf("\n");
            flint_printf("a (= b^n)    = "), fmpz_print(a), flint_printf("\n");
            flint_printf("b            = "), fmpz_print(b), flint_printf("\n");
            flint_printf("c (= (a)^(1/n) = "), fmpz_print(c), flint_printf("\n");
            flint_printf("d (= c^n)    = "), fmpz_print(d), flint_printf("\n");
            flint_printf("n            = %wd\n", n);
            flint_printf("ans          = %d\n", ans);
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(d);
        fmpz_clear(p);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
