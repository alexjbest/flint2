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

    Copyright (C) 2016 Alex J. Best

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_sparse_mat.h"

void
nmod_sparse_mat_set(nmod_sparse_mat_t B, const nmod_sparse_mat_t A)
{
    slong i;

    if (B == A || A->c == 0)
        return;

    for (i = 0; i < A->r; i++)
    {
        B->row_alloc[i] = A->row_alloc[i];
        B->row_supports[i] = A->row_supports[i];
        B->rows[i] = flint_malloc(B->row_alloc[i] * sizeof(nmod_sparse_mat_entry_struct));
        flint_mpn_copyi(B->rows[i], A->rows[i], 2 * B->row_alloc[i]);
    }
}
