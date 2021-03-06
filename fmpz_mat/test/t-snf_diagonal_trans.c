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

    Copyright (C) 2015 Alex J. Best

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mat.h"

int
main(void)
{
    slong iter;
    FLINT_TEST_INIT(state);

    flint_printf("snf_diagonal_trans....");
    fflush(stdout);

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        fmpz_mat_t A, S, S2, U, V;
        fmpz_t det;
        slong b, i, m, n;
        int equal;

        m = n_randint(state, 10);
        n = n_randint(state, 10);
        b = 1 + n_randint(state, 10) * n_randint(state, 10);

        fmpz_mat_init(A, m, n);
        fmpz_mat_init(S, m, n);
        fmpz_mat_init(S2, m, n);
        fmpz_mat_init(U, m, m);
        fmpz_mat_init(V, n, n);

        for (i = 0; i < FLINT_MIN(m, n); i++)
        {
            fmpz_randtest_unsigned(fmpz_mat_entry(A, i, i), state, b);
        }

        fmpz_mat_snf_diagonal_trans(S, U, V, A);

        if (!fmpz_mat_is_in_snf(S))
        {
            flint_printf("FAIL:\n");
            flint_printf("matrix not in snf!\n");
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(S); flint_printf("\n\n");
            abort();
        }

        fmpz_init(det);

        fmpz_mat_det(det, U);
        if (!fmpz_is_pm1(det))
        {
            flint_printf("FAIL:\n");
            flint_printf("transformation matrices should have determinant +-1, U does not!\n");
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(U); flint_printf("\n\n");
            fmpz_mat_print_pretty(V); flint_printf("\n\n");
            fmpz_mat_print_pretty(S); flint_printf("\n\n");
            abort();
        }

        fmpz_mat_det(det, V);
        if (!fmpz_is_pm1(det))
        {
            flint_printf("FAIL:\n");
            flint_printf("transformation matrices should have determinant +-1, V does not!\n");
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(U); flint_printf("\n\n");
            fmpz_mat_print_pretty(V); flint_printf("\n\n");
            fmpz_mat_print_pretty(S); flint_printf("\n\n");
            abort();
        }

        fmpz_clear(det);

        fmpz_mat_mul(S2, U, A);
        fmpz_mat_mul(S2, S2, V);
        equal = fmpz_mat_equal(S, S2);

        if (!equal)
        {
            flint_printf("FAIL:\n");
            flint_printf("multiplying by the transformation matrix should give the same SNF!\n");
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(U); flint_printf("\n\n");
            fmpz_mat_print_pretty(V); flint_printf("\n\n");
            fmpz_mat_print_pretty(S); flint_printf("\n\n");
            fmpz_mat_print_pretty(S2); flint_printf("\n\n");
            abort();
        }

        fmpz_mat_snf(S2, A);
        equal = fmpz_mat_equal(S, S2);

        if (!equal)
        {
            flint_printf("FAIL:\n");
            flint_printf("snfs produced by different methods should be the same!\n");
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(S); flint_printf("\n\n");
            fmpz_mat_print_pretty(S2); flint_printf("\n\n");
            abort();
        }
        fmpz_mat_clear(U);
        fmpz_mat_clear(V);
        fmpz_mat_clear(S2);
        fmpz_mat_clear(S);
        fmpz_mat_clear(A);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}

