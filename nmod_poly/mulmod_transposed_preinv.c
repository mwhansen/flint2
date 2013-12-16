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
_nmod_poly_mulmod_transposed_preinv(mp_ptr rop,
                                    mp_srcptr b, slong lenb,
                                    mp_srcptr a,
                                    mp_srcptr f, slong lenf,
                                    mp_srcptr finv, slong lenfinv,
                                    const nmod_t p)
{
    mp_limb_t *temp, *t1, *t2, *t3;

    temp = _nmod_vec_init(lenf);
    t1 = _nmod_vec_init(2 * lenf - 2);
    t2 = _nmod_vec_init(2 * lenf - 2);
    t3 = _nmod_vec_init(lenf - 1);
    
    /* t1 = rev(f) * a div x^n */
    _nmod_poly_reverse(temp, f, lenf, lenf);
    _nmod_poly_mul(t1, temp, lenf, a, lenf - 1, p);
    _nmod_poly_shift_right(t1, t1, 2*lenf - 2, lenf - 1);

    /* (t2 = ) rop = rev(b) * a div x^(n-1) */
    _nmod_poly_reverse(temp, b, lenb, lenf - 1);
    _nmod_poly_mul(t2, temp, lenf - 1, a, lenf - 1, p);
    _nmod_poly_shift_right(rop, t2, 2*lenf - 3, lenf - 2);

    /* t3 = x * (rev(f)^(-1) * rev(b) * t1 mod x^(n-1)) */
    _nmod_poly_mullow(t2, temp, lenf - 1, finv, lenfinv, lenf - 2, p);
    _nmod_poly_mullow(t3 + 1, t2, lenf - 2, t1, lenf - 1, lenf - 2, p);

    /* rop = rop - x*t3 */
    _nmod_poly_sub(rop, rop, lenf - 1, t3, lenf - 1, p);

    _nmod_vec_clear(temp);
    _nmod_vec_clear(t1);
    _nmod_vec_clear(t2);
    _nmod_vec_clear(t3);
}

/* Assume rop and a have room for f->length - 1 coeffs */
void
nmod_poly_mulmod_transposed_preinv(mp_ptr rop,
                                   const nmod_poly_t b,
                                   mp_srcptr a,
                                   const nmod_poly_t f,
                                   const nmod_poly_t finv)
{
    
    _nmod_poly_mulmod_transposed_preinv(rop,
                                        b->coeffs, b->length,
                                        a,
                                        f->coeffs, f->length,
                                        finv->coeffs, finv->length,
                                        f->mod);
}
