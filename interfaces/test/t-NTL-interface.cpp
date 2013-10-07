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

    Copyright (C) 2011 William Hart

******************************************************************************/

#include <NTL/ZZ.h>
#include <NTL/ZZX.h>

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_mod_poly.h"
#include "ulong_extras.h"
#include "NTL-interface.h"

NTL_CLIENT

int test_ZZ_to_fmpz()
{
    int i, result;
    flint_rand_t state;
    mp_bitcnt_t bits, randbits;
    fmpz_t int1, int2;
   
    ZZ z;

    printf("ZZ_to_fmpz....");
    fflush(stdout);

    flint_randinit(state);

    /* Check conversion */
    for (i = 0; i < 10000; i++)
    {
        bits = n_randint(state, 1000) + 1;
        randbits = n_randint(state, bits);
      
        fmpz_init(int1);
        fmpz_init(int2);
        
        fmpz_randbits(int1, state, randbits);

        fmpz_get_ZZ(z, int1);
        fmpz_set_ZZ(int2, z);
      
        result = fmpz_equal(int1, int2);
        if (!result)
        {
           printf("FAIL:\n");
           printf("int1 = %ld  ", *int1); fmpz_print(int1); printf("\n");
           printf("int2 = %ld  ", *int2); fmpz_print(int2); printf("\n");
           return 0;
        }

        fmpz_clear(int1);
        fmpz_clear(int2);
    }    

    return 1;
}

int test_ZZX_to_fmpz_poly()
{
    fmpz_poly_t f_poly1, f_poly2;
    slong length;
    mp_bitcnt_t bits;
    flint_rand_t state;
    int i, result;
   
    printf("ZZX_to_fmpz_poly....");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing of a and c */
    for (i = 0; i < 10000; i++)
    {
        bits = n_randint(state, 1000) + 1;
        length = n_randint(state, 1000);

        fmpz_poly_init(f_poly1);
        fmpz_poly_init(f_poly2);
        
        ZZX ZZX_poly;
          
        fmpz_poly_randtest(f_poly1, state, length, bits);
          
        fmpz_poly_get_ZZX(ZZX_poly, f_poly1);
        fmpz_poly_set_ZZX(f_poly2, ZZX_poly);
          
        result = fmpz_poly_equal(f_poly1, f_poly2);  
        if (!result)
        {
           printf("FAIL:\n");
           printf("f_poly1 = "); fmpz_poly_print(f_poly1); printf("\n");
           printf("f_poly2 = "); fmpz_poly_print(f_poly2); printf("\n");
           return 0;
        }
          
        fmpz_poly_clear(f_poly1);
        fmpz_poly_clear(f_poly2);
    }
      
    return 1;
}

int test_ZZ_pX_to_fmpz_mod_poly()
{
    fmpz_t p;
    fmpz_mod_poly_t f_poly1, f_poly2;
    slong length;
    flint_rand_t state;
    int i, result;
   
    printf("ZZ_pX_to_fmpz_mod_poly....");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing of a and c */
    for (i = 0; i < 10000; i++)
    {
        ZZ_pX ZZ_pX_poly;
        ZZ mod;

        length = n_randint(state, 1000);

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);

        fmpz_get_ZZ(mod, p);
        ZZ_p::init(mod);

        fmpz_mod_poly_init(f_poly1, p);
        fmpz_mod_poly_init(f_poly2, p);
          
        fmpz_mod_poly_randtest(f_poly1, state, length);
          
        fmpz_mod_poly_get_ZZ_pX(ZZ_pX_poly, f_poly1);
        fmpz_mod_poly_set_ZZ_pX(f_poly2, ZZ_pX_poly);
          
        result = fmpz_mod_poly_equal(f_poly1, f_poly2);  
        if (!result)
        {
           printf("FAIL:\n");
           printf("f_poly1 = "); fmpz_mod_poly_print(f_poly1); printf("\n");
           printf("f_poly2 = "); fmpz_mod_poly_print(f_poly2); printf("\n");
           return 0;
        }
        
        fmpz_clear(p);
        fmpz_mod_poly_clear(f_poly1);
        fmpz_mod_poly_clear(f_poly2);
    }
      
    return 1;
}

int test_ZZ_pE_to_fq()
{
    fmpz_t p;
    fq_t f1, f2;
    flint_rand_t state;
    int i, result;

    fq_ctx_t ctx;
   
    printf("ZZ_pE_to_fq....");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing of a and c */
    for (i = 0; i < 10000; i++)
    {
        fmpz_mod_poly_t fmod;
        long d;
        ZZ prime;
        ZZ_pX mod;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 6), 1));

        d = n_randint(state, 10) + 1;

        fmpz_get_ZZ(prime, p);
        ZZ_p::init(prime);

        BuildIrred(mod, d);
        ZZ_pE::init(mod);

        fmpz_mod_poly_init(fmod, p);
        fmpz_mod_poly_set_ZZ_pX(fmod, mod);

        fq_ctx_init_modulus(ctx, p, d, fmod, "a");

        ZZ_pE zzpe;

        fq_init(f1);
        fq_init(f2);

        fq_randtest(f1, state, ctx);

        fq_get_ZZ_pE(zzpe, f1);
        fq_set_ZZ_pE(f2, zzpe);
          
        result = fq_equal(f1, f2);
        if (!result)
        {
           printf("FAIL:\n");
           printf("p = "); fmpz_print(p); printf("\n");
           printf("mod = "); fmpz_mod_poly_print_pretty(fmod, "x"); printf("\n");
           printf("f1 = "); fq_print_pretty(f1, ctx); printf(" - %ld", f1->length); printf("\n");
           printf("zzpe:"); cout << zzpe; printf("\n");
           printf("f2 = "); fq_print_pretty(f2, ctx); printf(" - %ld", f2->length); printf("\n");
           return 0;
        }
        
        fmpz_clear(p);
        fq_clear(f1);
        fq_clear(f2);

        fmpz_mod_poly_clear(fmod);
        fq_ctx_clear(ctx);
    }
      
    return 1;
}

int test_ZZ_pEX_to_fq_poly()
{
    fmpz_t p;
    fq_poly_t f1, f2;
    slong length;
    flint_rand_t state;
    int i, result;

    fq_ctx_t ctx;
   
    printf("ZZ_pEX_to_fq_poly....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 10000; i++)
    {
        fmpz_mod_poly_t fmod;
        long d;
        ZZ prime;
        ZZ_pX mod;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 6), 1));

        d = n_randint(state, 10) + 1;

        fmpz_get_ZZ(prime, p);
        ZZ_p::init(prime);

        BuildIrred(mod, d);
        ZZ_pE::init(mod);

        fmpz_mod_poly_init(fmod, p);
        fmpz_mod_poly_set_ZZ_pX(fmod, mod);

        fq_ctx_init_modulus(ctx, p, d, fmod, "a");

        ZZ_pEX zzpex;

        fq_poly_init(f1);
        fq_poly_init(f2);

        length = n_randint(state, 1000);

        fq_poly_randtest(f1, state, length, ctx);

        fq_poly_get_ZZ_pEX(zzpex, f1);
        fq_poly_set_ZZ_pEX(f2, zzpex);
          
        result = fq_poly_equal(f1, f2);
        if (!result)
        {
           printf("FAIL:\n");
           printf("p = "); fmpz_print(p); printf("\n");
           printf("mod = "); fmpz_mod_poly_print_pretty(fmod, "x"); printf("\n");
           printf("f1 = "); fq_poly_print_pretty(f1, "x", ctx); printf("\n");
           printf("zzpex:"); cout << zzpex; printf("\n");
           printf("f2 = "); fq_poly_print_pretty(f2, "x", ctx); printf("\n");
           return 0;
        }
        
        fmpz_clear(p);
        fq_poly_clear(f1);
        fq_poly_clear(f2);

        fmpz_mod_poly_clear(fmod);
        fq_ctx_clear(ctx);
    }
      
    return 1;
}


int
main(void)
{
    int r = 1;
    
    if ((r &= test_ZZ_to_fmpz())) printf("PASS\n");
    if ((r &= test_ZZX_to_fmpz_poly())) printf("PASS\n");
    if ((r &= test_ZZ_pX_to_fmpz_mod_poly())) printf("PASS\n");
    if ((r &= test_ZZ_pE_to_fq())) printf("PASS\n");
    if ((r &= test_ZZ_pEX_to_fq_poly())) printf("PASS\n");

    if (!r) abort();

    _fmpz_cleanup();
    
    return 0;
}
