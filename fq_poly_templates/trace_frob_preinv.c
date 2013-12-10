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
TEMPLATE(T, poly_trace_frob_low_preinv)(TEMPLATE(T, poly_t) rop,
                                        const TEMPLATE(T, poly_t) alpha,
                                        slong m,
                                        const TEMPLATE(T, poly_t) beta,
                                        const TEMPLATE(T, poly_t) f,
                                        const TEMPLATE(T, poly_t) finv,
                                        const TEMPLATE(T, ctx_t) ctx)
{
    int i, l = (int)(FLINT_BIT_COUNT(m) - 1);
    TEMPLATE(T, poly_t) tau, tau_prime, mu, mu_prime, temp;
    

    TEMPLATE(T, poly_init)(tau, ctx);
    TEMPLATE(T, poly_init)(tau_prime, ctx);
    TEMPLATE(T, poly_init)(mu, ctx);
    TEMPLATE(T, poly_init)(mu_prime, ctx);
    TEMPLATE(T, poly_init)(temp, ctx);
                                        
    /* Stage 0 */
    TEMPLATE(T, poly_compose_mod_preinv)(tau, alpha, beta, f, finv, ctx);
    TEMPLATE(T, poly_set)(mu, beta, ctx);
    if (m & WORD(1))
    {
        TEMPLATE(T, poly_add)(tau_prime, alpha, tau, ctx);
        TEMPLATE(T, poly_set)(mu_prime, beta, ctx);
    }
    else
    {
        TEMPLATE(T, poly_set)(tau_prime, alpha, ctx);
        TEMPLATE(T, poly_gen)(mu_prime, ctx);
    }
    for (i = 1; i < l + 1; i++)
    {
        TEMPLATE(T, poly_compose_mod_preinv)(temp, tau, mu, f, finv, ctx);
        TEMPLATE(T, poly_add)(tau, tau, temp, ctx);
        if (m & WORD(1) << i)
        {
            TEMPLATE(T, poly_compose_mod_preinv)(temp, tau, mu_prime, f, finv, ctx);
            TEMPLATE(T, poly_add)(tau_prime, tau_prime, temp, ctx);
        }

        TEMPLATE(T, poly_compose_mod_preinv)(temp, mu, mu, f, finv, ctx);
        TEMPLATE(T, poly_set)(mu, temp, ctx);
        if (m & WORD(1) << i)
        {
            TEMPLATE(T, poly_compose_mod_preinv)(temp, mu, mu_prime, f, finv, ctx);
            TEMPLATE(T, poly_set)(mu_prime, temp, ctx);
        }
    }

    TEMPLATE(T, poly_set)(rop, tau_prime, ctx);

    TEMPLATE(T, poly_clear)(tau, ctx);
    TEMPLATE(T, poly_clear)(tau_prime, ctx);
    TEMPLATE(T, poly_clear)(mu, ctx);
    TEMPLATE(T, poly_clear)(mu_prime, ctx);
    TEMPLATE(T, poly_clear)(temp, ctx);
}


void
TEMPLATE(T, poly_trace_frob_preinv)(TEMPLATE(T, poly_t) rop,
                                    const TEMPLATE(T, poly_t) alpha,
                                    slong m,
                                    const TEMPLATE(T, poly_t) beta,
                                    const TEMPLATE(T, poly_t) f,
                                    const TEMPLATE(T, poly_t) finv,
                                    const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, poly_trace_frob_low_preinv)(rop, alpha, m, beta, f, finv, ctx);
}



#endif
