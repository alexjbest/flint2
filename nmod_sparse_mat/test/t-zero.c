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

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_sparse_mat.h"
#include "ulong_extras.h"

int
main(void)
{
    slong m, n, i, j, rep, N;
    FLINT_TEST_INIT(state);

    flint_printf("zero....");
    fflush(stdout);

    

    for (rep = 0; rep < 100 * flint_test_multiplier(); rep++)
    {
        nmod_sparse_mat_t A;

        m = n_randint(state, 20);
        n = n_randint(state, 20);
        N = n_randint(state, 20);

        nmod_sparse_mat_init(A, m, n, N);

	/*TODO randtest here */
        nmod_sparse_mat_zero(A);

        for (i = 0; i < m; i++)
        {
            for (j = 0; j < n; j++)
            {
                if (nmod_sparse_mat_get_entry(A, i, j) != 0)
                {
                    flint_printf("FAIL: nonzero entry\n");
                    abort();
                }
            }
        }

        nmod_sparse_mat_clear(A);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
