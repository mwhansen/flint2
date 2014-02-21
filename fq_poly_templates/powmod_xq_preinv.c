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
__TEMPLATE(T, poly_powmod_xq_preinv_compose)(TEMPLATE(T, poly_t) rop,
                                             const TEMPLATE(T, poly_t) f,
                                             const TEMPLATE(T, poly_t) finv,
                                             const TEMPLATE(T, ctx_t) ctx)
{
    int j, bc;
    slong ell, d = TEMPLATE(T, ctx_degree)(ctx);
    TEMPLATE(T, t) zeta_1, zeta;
    TEMPLATE(T, poly_t) xi_1, temp;
    TEMPLATE(T, poly_struct) *xi;

    /* xi_1 = x ^ p mod f */
    xi = rop;
    TEMPLATE(T, poly_init)(xi_1, ctx);
    TEMPLATE(T, poly_gen)(xi, ctx);
    TEMPLATE(T, poly_powmod_fmpz_sliding_preinv)(xi_1, xi,
                                                 TEMPLATE(T, ctx_prime)(ctx),
                                                 0, f, finv, ctx);
    TEMPLATE(T, poly_set)(xi, xi_1, ctx);

    /* zeta_1 = zeta = a^p */
    TEMPLATE(T, init)(zeta, ctx);
    TEMPLATE(T, init)(zeta_1, ctx);
    TEMPLATE(T, gen)(zeta, ctx);
    TEMPLATE(T, pow)(zeta_1, zeta, TEMPLATE(T, ctx_prime)(ctx), ctx);
    TEMPLATE(T, set)(zeta, zeta_1, ctx);

    TEMPLATE(T, poly_init)(temp, ctx);
    bc = (int)FLINT_BIT_COUNT(d) - 2;
    for (j = bc; j >= 0; j--)
    {
        TEMPLATE(T, poly_fit_length)(temp, xi->length, ctx);
        for (ell = 0; ell < xi->length; ell++)
        {
            TEMPLATE(T, compose)(temp->coeffs + ell,
                                 xi->coeffs + ell,
                                 zeta, ctx);
        }
        _TEMPLATE(T, poly_set_length)(temp, xi->length, ctx);
        _TEMPLATE(T, poly_normalise)(temp, ctx);
        TEMPLATE(T, poly_compose_mod_preinv)(xi, temp, xi, f, finv, ctx);

        TEMPLATE(T, compose)(zeta, zeta, zeta, ctx);

        if (d & (WORD(1) << j))
        {
            TEMPLATE(T, poly_fit_length)(temp, xi->length, ctx);
            for (ell = 0; ell < xi->length; ell++)
            {
                TEMPLATE(T, compose)(temp->coeffs + ell,
                                     xi->coeffs + ell,
                                     zeta_1, ctx);
            }
            _TEMPLATE(T, poly_set_length)(temp, xi->length, ctx);
            _TEMPLATE(T, poly_normalise)(temp, ctx);
            TEMPLATE(T, poly_compose_mod_preinv)(xi, temp, xi_1, f, finv, ctx);

            TEMPLATE(T, compose)(zeta, zeta, zeta_1, ctx);
        }
    }

    TEMPLATE(T, poly_clear)(temp, ctx);
    TEMPLATE(T, poly_clear)(xi_1, ctx);
    TEMPLATE(T, clear)(zeta_1, ctx);
    TEMPLATE(T, clear)(zeta, ctx);
}


void
__TEMPLATE(T, poly_powmod_xq_preinv_direct)(TEMPLATE(T, poly_t) rop,
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

void
TEMPLATE(T, poly_powmod_xq_preinv)(TEMPLATE(T, poly_t) rop,
                                   const TEMPLATE(T, poly_t) f,
                                   const TEMPLATE(T, poly_t) finv,
                                   const TEMPLATE(T, ctx_t) ctx)
{
    if (TEMPLATE(CAP_T, POLY_POWMOD_XQ_PREINV_USE_DIRECT)(f, ctx))
        __TEMPLATE(T, poly_powmod_xq_preinv_direct)(rop, f, finv, ctx);
    else
        __TEMPLATE(T, poly_powmod_xq_preinv_compose)(rop, f, finv, ctx);
}

#endif
