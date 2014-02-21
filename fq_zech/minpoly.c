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

#include "fq_zech.h"

void
fq_zech_minpoly(nmod_poly_t g, const fq_zech_t op, const fq_zech_ctx_t ctx)
{
    fq_nmod_t nop;
    fq_nmod_init(nop, ctx->fq_nmod_ctx);
    fq_zech_get_fq_nmod(nop, op, ctx);
    fq_nmod_minpoly(g, nop, ctx->fq_nmod_ctx);
    fq_nmod_clear(nop, ctx->fq_nmod_ctx);
}
