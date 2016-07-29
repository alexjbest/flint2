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
nmod_sparse_mat_init_set(nmod_sparse_mat_t mat, const nmod_sparse_mat_t src)
{
    slong rows = src->r;
    slong cols = src->c;

    if ((rows) && (cols))
    {
        slong i;
        mat->rows = flint_malloc(rows * sizeof(nmod_sparse_mat_entry_struct *));
        mat->row_alloc = flint_malloc(rows * sizeof(mp_limb_t));
        mat->row_supports = flint_malloc(rows * sizeof(mp_limb_t));

        for (i = 0; i < rows; i++)
        {
            mat->row_alloc[i] = src->row_alloc[i];
            mat->row_supports[i] = src->row_supports[i];
            mat->rows[i] = flint_malloc(mat->row_alloc[i] * sizeof(nmod_sparse_mat_entry_struct));
            flint_mpn_copyi(mat->rows[i], src->rows[i], 2 * mat->row_alloc[i]);
        }
    }
    else
        mat->rows = NULL;

    mat->r = rows;
    mat->c = cols;

    mat->mod = src->mod;
}
