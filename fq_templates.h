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
TEMPLATE(T, ctx_init_minpoly)(TEMPLATE(T, ctx_t) ctxp, const TEMPLATE(T, t) op,
                                      const char *var, const TEMPLATE(T, ctx_t) ctx);

void
TEMPLATE(T, gcdinv)(TEMPLATE(T, t) rop, TEMPLATE(T, t) inv,
                    const TEMPLATE(T, t) op,
                    const TEMPLATE(T, ctx_t) ctx);

int
TEMPLATE(T, is_invertible)(const TEMPLATE(T, t) op,
                           const TEMPLATE(T, ctx_t) ctx);

int
TEMPLATE(T, is_invertible_f)(TEMPLATE(T, t) rop, const TEMPLATE(T, t) op,
                             const TEMPLATE(T, ctx_t) ctx);

void
TEMPLATE(T, div)(TEMPLATE(T, t) rop, const TEMPLATE(T, t) op1,
                 const TEMPLATE(T, t) op2, const TEMPLATE(T, ctx_t) ctx);

/* Root finding */

void
TEMPLATE(T, compose)(TEMPLATE(T, t) rop,
                     const TEMPLATE(T, t) op, const TEMPLATE(T, t) xp,
                     const TEMPLATE(T, ctx_t) ctx);

int
TEMPLATE(T, root_edf)(TEMPLATE(T, t) rop, const TEMPLATE(T, t) op, slong n,
                      const TEMPLATE(T, ctx_t) ctx);

int
TEMPLATE(T, root_ds)(TEMPLATE(T, t) rop, const TEMPLATE(T, t) op, slong n,
                     const TEMPLATE(T, ctx_t) ctx);

void
__TEMPLATE(T, root_ds_s_ell)(fmpz_t s, fmpz_t ell, fmpz_t e, slong t,
                             const TEMPLATE(T, ctx_t) ctx);

void
__TEMPLATE(T, root_ds_qm1ot)(TEMPLATE(T, t) rop, const TEMPLATE(T, t) op,
                             const fmpz_t s, const fmpz_t ell, const fmpz_t e,
                             slong t, const TEMPLATE(T, ctx_t) ctx);

void
__TEMPLATE(T, root_ds_xi1zeta1)(TEMPLATE(T, t) xi1, TEMPLATE(T, t) zeta1,
                                TEMPLATE(T, t) xp,
                                const TEMPLATE(T, t) op,
                                const fmpz_t s,
                                const TEMPLATE(T, ctx_t) ctx);

void
__TEMPLATE(T, root_ds_xizeta)(TEMPLATE(T, t) xi, TEMPLATE(T, t) zeta,
                              ulong i,
                              const TEMPLATE(T, t) xi1,
                              const TEMPLATE(T, t) zeta1,
                              const TEMPLATE(T, ctx_t) ctx);

int
_TEMPLATE(T, root_prime)(TEMPLATE(T, t) rop, const TEMPLATE(T, t) op, slong n,
                         const TEMPLATE(T, ctx_t) ctx);

#endif
