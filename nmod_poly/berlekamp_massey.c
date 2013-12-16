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

void
nmod_poly_berlekamp_massey(nmod_poly_t rop, mp_srcptr c, slong k)
{
    slong i, j;
    slong L = 0, m = 0;
    nmod_poly_t lambda_rev, B, Bshift;
    mp_limb_t delta_i, delta, temp;

    nmod_poly_init2_preinv(lambda_rev, rop->mod.n, rop->mod.ninv, k + 1);
    nmod_poly_init2_preinv(B, rop->mod.n, rop->mod.ninv, k + 1);
    nmod_poly_init2_preinv(Bshift, rop->mod.n, rop->mod.ninv, k + 1);

    nmod_poly_set_coeff_ui(lambda_rev, 0, 1);

    delta = 1;
    
    for (i = 1; i <= 2*k; i++)
    {
        delta_i = 0;
        for (j = 0; j < lambda_rev->length; j++)
        {
            delta_i = n_addmod(delta_i,
                               n_mulmod2_preinv(lambda_rev->coeffs[j],
                                                c[i - j - 1],
                                                rop->mod.n, rop->mod.ninv),
                               rop->mod.n);
        }
         
        if (delta_i == 0)
        {
            m += 1;
        }
        else if (2 * L < i)
        {
            /* temp = delta_i / delta */
            temp = n_invmod(delta, rop->mod.n);
            temp = n_mulmod2_preinv(temp, delta_i, rop->mod.n, rop->mod.ninv);
            
            nmod_poly_shift_left(Bshift, B, m + 1);
            nmod_poly_scalar_mul_nmod(Bshift, Bshift, temp);

            nmod_poly_set(B, lambda_rev);

            nmod_poly_sub(lambda_rev, lambda_rev, Bshift);
            
            m = 0;
            L = i - L;
            delta = delta_i;
        }
        else
        {
            /* temp = delta_i / delta */
            temp = n_invmod(delta, rop->mod.n);
            temp = n_mulmod2_preinv(temp, delta_i, rop->mod.n, rop->mod.ninv);

            nmod_poly_shift_left(Bshift, B, m + 1);
            nmod_poly_scalar_mul_nmod(Bshift, Bshift, temp);

            nmod_poly_sub(lambda_rev, lambda_rev, Bshift);

            m += 1;
        }
    }

    nmod_poly_reverse(rop, lambda_rev, L + 1);

    nmod_poly_clear(lambda_rev);
    nmod_poly_clear(B);
    nmod_poly_clear(Bshift);
}
