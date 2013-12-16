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

void
fmpz_mod_poly_berlekamp_massey(fmpz_mod_poly_t rop, const fmpz *c, ulong k)
{
    ulong i, j;
    ulong L = 0, m = 0;
    fmpz_mod_poly_t lambda_rev, B, Bshift;
    fmpz_t delta_i, delta, temp;

    fmpz_mod_poly_init2(lambda_rev, &(rop->p), k + 1);
    fmpz_mod_poly_init2(B, &(rop->p), k + 1);
    fmpz_mod_poly_init2(Bshift, &(rop->p), k + 1);

    fmpz_mod_poly_set_coeff_ui(lambda_rev, 0, 1);

    fmpz_init_set_ui(delta, 1);
    fmpz_init(delta_i);
    fmpz_init(temp);
    
    for (i = 1; i <= 2*k; i++)
    {
        fmpz_zero(delta_i);
        for (j = 0; j < lambda_rev->length; j++)
        {
            fmpz_addmul(delta_i, lambda_rev->coeffs + j, c + (i - j - 1));
        }
        fmpz_mod(delta_i, delta_i, &(rop->p));
         
        if (fmpz_is_zero(delta_i))
        {
            m += 1;
        }
        else if (2 * L < i)
        {
            /* temp = delta_i / delta */
            fmpz_invmod(temp, delta, &(rop->p));
            fmpz_mul(temp, delta_i, temp);
            
            fmpz_mod_poly_shift_left(Bshift, B, m + 1);
            fmpz_mod_poly_scalar_mul_fmpz(Bshift, Bshift, temp);

            fmpz_mod_poly_set(B, lambda_rev);

            fmpz_mod_poly_sub(lambda_rev, lambda_rev, Bshift);
            
            m = 0;
            L = i - L;
            fmpz_set(delta, delta_i);
        }
        else
        {
            /* temp = delta_i / delta */
            fmpz_invmod(temp, delta, &(rop->p));
            fmpz_mul(temp, delta_i, temp);
            fmpz_mod(temp, temp, &(rop->p));

            fmpz_mod_poly_shift_left(Bshift, B, m + 1);
            fmpz_mod_poly_scalar_mul_fmpz(Bshift, Bshift, temp);

            fmpz_mod_poly_sub(lambda_rev, lambda_rev, Bshift);

            m += 1;
        }
    }

    fmpz_mod_poly_reverse(rop, lambda_rev, L + 1);

    fmpz_mod_poly_clear(lambda_rev);
    fmpz_mod_poly_clear(B);
    fmpz_mod_poly_clear(Bshift);

    fmpz_clear(delta);
    fmpz_clear(delta_i);
    fmpz_clear(temp);
}
