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
_fmpz_mod_poly_mulmod_transposed_preinv(fmpz *rop,
                                        const fmpz *b, slong lenb,
                                        const fmpz *a,
                                        const fmpz *f, slong lenf,
                                        const fmpz *finv, slong lenfinv,
                                        const fmpz_t p)
{
    fmpz *temp, *t1, *t2, *t3;

    temp = _fmpz_vec_init(lenf);
    t1 = _fmpz_vec_init(2 * lenf - 2);
    t2 = _fmpz_vec_init(2 * lenf - 2);
    t3 = _fmpz_vec_init(lenf - 1);
    
    /* t1 = rev(f) * a div x^n */
    _fmpz_mod_poly_reverse(temp, f, lenf, lenf);
    _fmpz_mod_poly_mul(t1, temp, lenf, a, lenf - 1, p);
    _fmpz_mod_poly_shift_right(t1, t1, 2*lenf - 2, lenf - 1);

    /* (t2 = ) rop = rev(b) * a div x^(n-1) */
    _fmpz_mod_poly_reverse(temp, b, lenb, lenf - 1);
    _fmpz_mod_poly_mul(t2, temp, lenf - 1, a, lenf - 1, p);
    _fmpz_mod_poly_shift_right(rop, t2, 2*lenf - 3, lenf - 2);

    /* t3 = x * (rev(f)^(-1) * rev(b) * t1 mod x^(n-1)) */
    _fmpz_mod_poly_mullow(t2, temp, lenf - 1, finv, lenfinv, p, lenf - 2);
    _fmpz_mod_poly_mullow(t3 + 1, t2, lenf - 2, t1, lenf - 1, p, lenf - 2);

    /* rop = rop - x*t3 */
    _fmpz_mod_poly_sub(rop, rop, lenf - 1, t3, lenf - 1, p);

    _fmpz_vec_clear(temp, lenf);
    _fmpz_vec_clear(t1, 2 * lenf - 2);
    _fmpz_vec_clear(t2, 2 * lenf - 2);
    _fmpz_vec_clear(t3, lenf - 1);
}

/* Assume rop and a have room for f->length - 1 coeffs */
void
fmpz_mod_poly_mulmod_transposed_preinv(fmpz *rop,
                                       const fmpz_mod_poly_t b,
                                       const fmpz *a,
                                       const fmpz_mod_poly_t f,
                                       const fmpz_mod_poly_t finv)
{
    
    _fmpz_mod_poly_mulmod_transposed_preinv(rop,
                                            b->coeffs, b->length,
                                            a,
                                            f->coeffs, f->length,
                                            finv->coeffs, finv->length,
                                            &(f->p));
}
