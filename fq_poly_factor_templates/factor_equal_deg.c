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

    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010 William Hart
    Copyright (C) 2008 Richard Howell-Peak
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2013 Mike Hansen

******************************************************************************/


#ifdef T

#include "templates.h"

#include "ulong_extras.h"

void
TEMPLATE(T, poly_factor_equal_deg) (TEMPLATE(T, poly_factor_t) factors,
                                    const TEMPLATE(T, poly_t) pol, slong d,
                                    const TEMPLATE(T, ctx_t) ctx)
{
    if (pol->length == d + 1)
    {
        TEMPLATE(T, poly_factor_insert) (factors, pol, 1, ctx);
    }
    else
    {
        TEMPLATE(T, poly_t) f, g, h, r;
        int factors_found;
        flint_rand_t state;

        TEMPLATE(T, poly_init) (f, ctx);
        TEMPLATE(T, poly_init) (g, ctx);
        TEMPLATE(T, poly_init) (h, ctx);
        TEMPLATE(T, poly_init) (r, ctx);

        flint_randinit(state);

        factors_found = 0;
        while (!factors_found)
        {
            factors_found = TEMPLATE(T, poly_factor_equal_deg_prob)
                                (f, g, state, pol, d, ctx);
        }

        flint_randclear(state);

        TEMPLATE(T, poly_divrem) (h, r, pol, f, ctx);
        TEMPLATE(T, poly_factor_equal_deg) (factors, f, d, ctx);
        if (factors_found == 2)
        {
            TEMPLATE(T, poly_divrem) (h, r, h, g, ctx);
            TEMPLATE(T, poly_factor_equal_deg) (factors, g, d, ctx);
        }
        if (!TEMPLATE(T, poly_is_one)(h, ctx))
            TEMPLATE(T, poly_factor_equal_deg) (factors, h, d, ctx);

        TEMPLATE(T, poly_clear) (f, ctx);
        TEMPLATE(T, poly_clear) (g, ctx);
        TEMPLATE(T, poly_clear) (h, ctx);
        TEMPLATE(T, poly_clear) (r, ctx);
    }
}


#endif
