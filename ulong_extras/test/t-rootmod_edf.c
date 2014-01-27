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

#include "ulong_extras.h"

int main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("rootmod_edf....");
    fflush(stdout);

    for (i = 0; i < 100 * flint_test_multiplier(); i++) /* Test random integers */
    {
        mp_limb_t a, b, p;
        slong n;
        
        n = n_randint(state, 32) + 2;
        
        p = n_randtest_prime(state, 0);
        a = n_randtest(state) % p;

        b = n_rootmod_edf(a, n, p);

        result = (b == 0 || n_powmod2(b, n, p) == a);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("p = %wu\n", p);
            flint_printf("n = %wu\n", n);
            flint_printf("a = %wu\n", a);
            flint_printf("b = %wu\n", b);
            abort();
        }
    }

    for (i = 0; i < 100 * flint_test_multiplier(); i++) /* Test random n^th powers */
    {
        mp_limb_t a, b, p;
        slong n;

        p = n_randtest_prime(state, 0);

        do 
            b = n_randtest(state) % p;
        while (b == 0);

        n = n_randint(state, 32) + 2;

        a = n_powmod2(b, n, p);
        b = n_rootmod_edf(a, n, p);

        result = (n_powmod2(b, n, p) == a);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("p = %wu\n", p);
            flint_printf("n = %wu\n", n);
            flint_printf("a = %wu\n", a);
            flint_printf("b = %wu\n", b);
            abort();
        }
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
