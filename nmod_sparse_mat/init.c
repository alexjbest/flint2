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
nmod_sparse_mat_init(nmod_sparse_mat_t mat, slong rows, slong cols, mp_limb_t n)
{
    if ((rows) && (cols))
    {
        slong i;
        mat->rows = flint_malloc(rows * sizeof(nmod_sparse_mat_entry_struct *));
        mat->row_alloc = flint_calloc(rows, sizeof(mp_limb_t));
        mat->row_supports = flint_calloc(rows, sizeof(mp_limb_t));

        for (i = 0; i < rows; i++)
        {
            mat->rows[i] = flint_malloc(0 * sizeof(nmod_sparse_mat_entry_struct));
        }
    }
    else
        mat->rows = NULL;

    mat->r = rows;
    mat->c = cols;

    _nmod_sparse_mat_set_mod(mat, n);
}
