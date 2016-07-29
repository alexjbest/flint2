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
#include <stdio.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_sparse_mat.h"
#include "ulong_extras.h"

void
nmod_sparse_mat_print_pretty(const nmod_sparse_mat_t mat)
{
    slong i, j;
    int width;
    char fmt[FLINT_BITS + 5];

    flint_printf("<%wd x %wd integer matrix mod %wu>\n", mat->r, mat->c, mat->mod.n);

    if (!(mat->c) || !(mat->r))
        return;

    width = n_sizeinbase(mat->mod.n, 10);

    flint_sprintf(fmt, "%%%dwu", width);

    for (i = 0; i < mat->r; i++)
    {
        flint_printf("[%d/%d:", mat->row_supports[i], mat->row_alloc[i]);

        for (j = 0; j < mat->row_alloc[i]; j++)
        {
            flint_printf(fmt, mat->rows[i][j].pos);
            flint_printf(" ");
            flint_printf(fmt, mat->rows[i][j].entry);
            if (j + 1 < mat->c)
                flint_printf(" ");
        }

        flint_printf("]\n");
    }
}

