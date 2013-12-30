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

    Copyright (C) 2007, 2008 William Hart
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

******************************************************************************/

#include "fmpz_mod_poly.h"

void
_fmpz_mod_poly_resultant_euclidean(fmpz_t rop, const fmpz* poly1, slong len1, 
                                   const fmpz* poly2, slong len2,
                                   const fmpz_t p)
{
    if (poly1 == poly2)
    {
        fmpz_zero(rop);
        return;
    }
    else if (len2 == 1)
    {
        if (len1 == 1)
        {
            fmpz_one(rop);
        }
        else if (len1 == 2)
        {
            fmpz_set(rop, poly2);
        }
        else
        {
            fmpz_powm_ui(rop, poly2, len1 - 1, p);
        }
        return;
    }
    else  /* len1 >= len2 >= 2 */
    {
        fmpz *u, *v, *r, *t, *w;
        slong l0, l1, l2;
        fmpz_t lc, inv;

        fmpz_one(rop);
        fmpz_init(lc);
        fmpz_init(inv);

        w = _fmpz_vec_init(3 * len1);
        u = w;
        v = w + len1;
        r = v + len1;

        _fmpz_vec_set(u, poly1, len1);
        _fmpz_vec_set(v, poly2, len2);
        l1 = len1;
        l2 = len2;

        do
        {
            l0 = l1;
            l1 = l2;
            fmpz_set(lc, v + l1 - 1);

            fmpz_invmod(inv, lc, p);
            _fmpz_mod_poly_rem(r, u, l0, v, l1, inv, p);
            l2 = l1 - 1;
            FMPZ_VEC_NORM(r, l2);
            {
                t = u;
                u = v;
                v = r;
                r = t;
            }

            if (l2 >= 1) 
            {
                fmpz_powm_ui(lc, lc, l0 - l2, p);
                fmpz_mul(rop, rop, lc);
                fmpz_mod(rop, rop, p);

                if (((l0 | l1) & 1) == 0)
                {
                    if (fmpz_cmp_ui(rop, 0) > 0)
                        fmpz_sub(rop, p, rop);
                }  
            }
            else 
            {
                if (l1 == 1)
                {
                    fmpz_powm_ui(lc, lc, l0 - 1, p);
                    fmpz_mul(rop, rop, lc);
                    fmpz_mod(rop, rop, p);
                }
                else
                {
                    fmpz_zero(rop);
                }
            }
        }
        while (l2 > 0);

        _fmpz_vec_clear(w, 3 * len1);
        fmpz_clear(inv);
        fmpz_clear(lc);
    }
}

void
fmpz_mod_poly_resultant_euclidean(fmpz_t rop, const fmpz_mod_poly_t f,
                                  const fmpz_mod_poly_t g)
{
    const slong len1 = f->length;
    const slong len2 = g->length;

    if (len1 == 0 || len2 == 0)
    {
        fmpz_zero(rop);
    }
    else
    {
        if (len1 >= len2)
        {
            _fmpz_mod_poly_resultant_euclidean(rop, f->coeffs, len1, 
                                               g->coeffs, len2, &(f->p));
        }
        else
        {
            _fmpz_mod_poly_resultant_euclidean(rop, g->coeffs, len2, 
                                               f->coeffs, len1, &(f->p));

            if (((len1 | len2) & WORD(1)) == WORD(0))
                if (fmpz_cmp_ui(rop, 0) > 0)
                    fmpz_sub(rop, &(f->p), rop);
        }
    }
}

