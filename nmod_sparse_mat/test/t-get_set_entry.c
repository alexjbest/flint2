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
#include "nmod_vec.h"
#include "nmod_sparse_mat.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, j, result, r, c;

    FLINT_TEST_INIT(state);

    flint_printf("get/set_entry....");
    fflush(stdout);

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_sparse_mat_t a;

        nmod_sparse_mat_init(a, n_randint(state, 10) + 1, n_randint(state, 10) + 1, n_randint(state, 10000));

        for (j = 0; j < 10000; j++)
        {
            mp_limb_t in, out;
            in = n_randlimb(state);
            r = (slong)n_randint(state, a->r);
            c = (slong)n_randint(state, a->c);
            nmod_sparse_mat_set_entry(a, r, c, in);
            out = nmod_sparse_mat_get_entry(a, r, c);

	    NMOD_RED(in, in, a->mod);
            result = (in == out);
            if (!result)
            {
                flint_printf("FAIL:\n\n");
		nmod_sparse_mat_print_pretty(a);
                flint_printf("r   = %wd\n\n", r);
                flint_printf("c   = %wd\n\n", c);
                flint_printf("in  = %wd\n\n", in);
                flint_printf("out = %wd\n\n", out);
                abort();
            }
        }

        nmod_sparse_mat_clear(a);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}
