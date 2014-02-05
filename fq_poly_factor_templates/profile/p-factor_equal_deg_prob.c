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



#include "flint.h"
#include <math.h>
#include "profiler.h"

#define nalgs 3
#define cpumin 2
#define ncases 10

int
main(int argc, char** argv)
{
    fmpz_t p;
    int c, n, reps = 0;
    slong d, lenf, m;
    TEMPLATE(T, ctx_t) ctx;
    TEMPLATE(T, poly_t) f;
    TEMPLATE(T, poly_factor_t) factors;
    slong i, l;
    double beta;

    double s[nalgs];
    
    FLINT_TEST_INIT(state);
    
    fmpz_init(p);
    fmpz_set_str(p, argv[1], 10);

    d = atol(argv[2]);
    lenf = atol(argv[3]);
    m = atol(argv[4]);

    TEMPLATE(T, ctx_init)(ctx, p, d, "a");

    TEMPLATE(T, poly_init)(f, ctx);
    
    for (c = 0; c < nalgs; c++)
        s[c] = 0.0;
            
    for (n = 0; n < ncases; n++)
    {
        double t[nalgs];
        int lo, loops = 1;

        /*
          Construct random elements of fq
        */
        {
            TEMPLATE(T, poly_randtest_irreducible)(f, state, lenf, ctx);
            TEMPLATE(T, poly_pow)(f, f, m, ctx);
        }
                
    loop:

        t[0] = 0.0;
        init_clock(0);
        prof_start();
        for (lo = 0; lo < loops; lo++)
        {
            TEMPLATE(T, poly_factor_init)(factors, ctx);
            TEMPLATE(T, poly_factor_equal_deg_cz)(factors, f, lenf - 1, ctx);
            TEMPLATE(T, poly_factor_clear)(factors, ctx);
        }
        prof_stop();
        t[0] += get_clock(0);

        
        t[1] = 0.0;
        init_clock(0);
        prof_start();
        for (lo = 0; lo < loops; lo++)
        {
            TEMPLATE(T, poly_factor_init)(factors, ctx);
            TEMPLATE(T, poly_factor_equal_deg_vzgs)(factors, f, lenf - 1, ctx);
            TEMPLATE(T, poly_factor_clear)(factors, ctx);
        }
        prof_stop();
        t[1] += get_clock(0);

        t[2] = 0.0;
        init_clock(0);
        prof_start();
        for (lo = 0; lo < loops; lo++)
        {
            TEMPLATE(T, poly_factor_init)(factors, ctx);
            TEMPLATE(T, poly_factor_equal_deg_ks)(factors, f, lenf - 1, ctx);
            TEMPLATE(T, poly_factor_clear)(factors, ctx);
        }
        prof_stop();
        t[2] += get_clock(0);

        for (c = 0; c < nalgs; c++)
            if (t[c] * FLINT_CLOCK_SCALE_FACTOR <= cpumin)
            {
                loops *= 10;
                goto loop;
            }
                
        for (c = 0; c < nalgs; c++)
            s[c] += t[c];
        reps += loops;
    }
            
    for (c = 0; c < nalgs; c++)
    {
        printf("%20f", s[c] / (double) reps);
        fflush(stdout);
    }
    printf("\n");
    
    TEMPLATE(T, poly_clear)(f, ctx);
    TEMPLATE(T, ctx_clear)(ctx);
    fmpz_clear(p);

    FLINT_TEST_CLEANUP(state);
    
    return 0;
}


#endif
