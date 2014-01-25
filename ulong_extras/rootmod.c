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

#include "nmod_poly.h"

mp_limb_t
n_rootmod(mp_limb_t a, slong n, mp_limb_t p)
{
    mp_limb_t found = 0;
    slong i;
    nmod_poly_struct *f;
    nmod_poly_t v, finv, x, xpmx, g;
    nmod_poly_factor_t sq_free, equal_deg;

    if (a == 0)
        return 0;

    nmod_poly_init(v, p);
    nmod_poly_init(x, p);
    nmod_poly_init(finv, p);
    nmod_poly_init(xpmx, p);
    nmod_poly_init(g, p);


    /* v = x^n - op */
    nmod_poly_fit_length(v, n + 1);
    _nmod_vec_zero(v->coeffs, n + 1);
    v->coeffs[n] = 1;
    v->coeffs[0] = n_negmod(a, p);
    _nmod_poly_set_length(v, n + 1);

    nmod_poly_gen(x);

    /* Squarefree factorization */
    nmod_poly_factor_init(sq_free);
    nmod_poly_factor_squarefree(sq_free, v);

    for (i = 0; i < sq_free->num; i++)
    {
        f = sq_free->p + i;
        nmod_poly_reverse(finv, f, f->length);
        nmod_poly_inv_series_newton(finv, finv, f->length);

        nmod_poly_powmod_x_ui_preinv(xpmx, p, f, finv);
        nmod_poly_sub(xpmx, xpmx, x);

        nmod_poly_gcd(g, xpmx, f);

        if (!nmod_poly_is_one(g))
        {
            nmod_poly_factor_init(equal_deg);
            nmod_poly_factor_equal_deg(equal_deg, g, 1);

            if (equal_deg->num > 0)
            {
                found = n_negmod(equal_deg->p->coeffs[0], p);
            }
            nmod_poly_factor_clear(equal_deg);

            if (found)
                break;
        }
    }

    nmod_poly_factor_clear(sq_free);
    nmod_poly_clear(v);
    nmod_poly_clear(x);
    nmod_poly_clear(finv);
    nmod_poly_clear(xpmx);
    nmod_poly_clear(g);

    return found;
}
