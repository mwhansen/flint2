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

#include "fmpz.h"
#include "fmpz_mod_poly.h"

int
fmpz_rootmod(fmpz_t rop, const fmpz_t op, slong n, const fmpz_t p)
{
    int found = 0;
    slong i;
    fmpz_mod_poly_struct *f;
    fmpz_mod_poly_t v, finv, x, xpmx, g;
    fmpz_mod_poly_factor_t sq_free, equal_deg;

    if (fmpz_is_zero(op))
    {
        fmpz_zero(rop);
        return 1;
    }

    fmpz_mod_poly_init(v, p);
    fmpz_mod_poly_init(x, p);
    fmpz_mod_poly_init(finv, p);
    fmpz_mod_poly_init(xpmx, p);
    fmpz_mod_poly_init(g, p);

    /* v = x^n - op */
    fmpz_mod_poly_fit_length(v, n + 1);
    fmpz_one(v->coeffs + n);
    fmpz_sub(v->coeffs, p, op);
    _fmpz_mod_poly_set_length(v, n + 1);

    fmpz_mod_poly_gen(x);

    /* Squarefree factorization */
    fmpz_mod_poly_factor_init(sq_free);
    fmpz_mod_poly_factor_squarefree(sq_free, v);

    for (i = 0; i < sq_free->num; i++)
    {
        f = sq_free->poly + i;
        fmpz_mod_poly_reverse(finv, f, f->length);
        fmpz_mod_poly_inv_series_newton(finv, finv, f->length);

        fmpz_mod_poly_powmod_x_fmpz_preinv(xpmx, p, f, finv);
        fmpz_mod_poly_sub(xpmx, xpmx, x);

        fmpz_mod_poly_gcd(g, xpmx, f);

        if (!fmpz_mod_poly_is_one(g))
        {
            fmpz_mod_poly_factor_init(equal_deg);
            fmpz_mod_poly_factor_equal_deg(equal_deg, g, 1);
            if (equal_deg->num > 0)
            {
                fmpz_sub(rop, p, equal_deg->poly->coeffs);
                found = 1;
            }
            fmpz_mod_poly_factor_clear(equal_deg);
            if (found)
                break;
        }
    }

    fmpz_mod_poly_factor_clear(sq_free);
    fmpz_mod_poly_clear(v);
    fmpz_mod_poly_clear(x);
    fmpz_mod_poly_clear(finv);
    fmpz_mod_poly_clear(xpmx);
    fmpz_mod_poly_clear(g);

    return found;
}
