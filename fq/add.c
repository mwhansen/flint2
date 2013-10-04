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

    Copyright (C) 2011, 2012 Sebastian Pancratz 
    Copyright (C) 2012 Andres Goens
 
******************************************************************************/

#include "fq.h"

void fq_add(fq_t rop, const fq_t op1, const fq_t op2, const fq_ctx_t ctx)
{
    long max = FLINT_MAX(op1->length, op2->length);

    fmpz_poly_fit_length(rop, max);

    _fmpz_mod_poly_add(rop->coeffs, 
                       op1->coeffs, op1->length, op2->coeffs, op2->length, 
                       fq_ctx_prime(ctx));

    _fmpz_poly_set_length(rop, max);
    _fmpz_poly_normalise(rop);
}
