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

    Copyright (C) 2012 Sebastian Pancratz

******************************************************************************/

#include "fq_poly.h"

void _fq_poly_derivative(fq_struct *rop, const fq_struct *op, long len, 
                                         const fq_ctx_t ctx)
{
    long i;

    for (i = 1; i < len; i++)
        fq_mul_ui(rop + (i - 1), op + i, i, ctx);
}

void fq_poly_derivative(fq_poly_t rop, const fq_poly_t op, const fq_ctx_t ctx)
{
    const long len = op->length;

    if (len < 2)
    {
        fq_poly_zero(rop);
    }
    else
    {
        fq_poly_fit_length(rop, len - 1);
        _fq_poly_derivative(rop->coeffs, op->coeffs, len, ctx);
        _fq_poly_set_length(rop, len - 1);
        _fq_poly_normalise(rop);
    }
}

