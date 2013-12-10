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


#ifdef T

#include "templates.h"

void
TEMPLATE(T, poly_powmod_xq_preinv)(TEMPLATE(T, poly_t) rop,
                                   const TEMPLATE(T, poly_t) f,
                                   const TEMPLATE(T, poly_t) finv,
                                   const TEMPLATE(T, ctx_t) ctx)
{
    fmpz_t q;
    TEMPLATE(T, poly_t) x;

    fmpz_init(q);
    TEMPLATE(T, ctx_order)(q, ctx);
    
    TEMPLATE(T, poly_init)(x, ctx);
    TEMPLATE(T, poly_gen)(x, ctx);
    TEMPLATE(T, poly_powmod_fmpz_sliding_preinv)(rop, x, q, 0, f, finv, ctx);

    TEMPLATE(T, poly_clear)(x, ctx);
    fmpz_clear(q);
}

#endif
