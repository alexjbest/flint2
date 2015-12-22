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

#include "fmpz_mat.h"

/* sets a to gcd(a,b) and b to lcm(a,b) using temporary fmpz_t t */
static void _gcdlcm(fmpz_t t, fmpz_t a, fmpz_t b)
{
    if (fmpz_equal(a, b)) return;
    fmpz_gcd(t, a, b);
    fmpz_divexact(b, b, t);
    fmpz_mul(b, b, a);
    fmpz_set(a, t);
}

void _fmpz_mat_snf_diagonal_trans(fmpz_mat_t S, fmpz_mat_t U, fmpz_mat_t V, const fmpz_mat_t A)
{
    fmpz_t d, x, y;
    slong i, j, k, n = FLINT_MIN(A->r, A->c);

    fmpz_init(d);
    fmpz_init(x);
    fmpz_init(y);
    fmpz_mat_set(S, A);
    for (i = 0; i < n; i++)
    {
        if (fmpz_sgn(fmpz_mat_entry(S, i, i)) < 0)
        {
            fmpz_neg(fmpz_mat_entry(S, i, i), fmpz_mat_entry(S, i, i));
            for (k = 0; k < A->r; k++)
                fmpz_neg(fmpz_mat_entry(U, i, k), fmpz_mat_entry(U, i, k));
        }
    }
    for (j = n - 1; j >= 0; j--)
    {
        for (i = 0; i < j; i++)
        {
            if (fmpz_equal(fmpz_mat_entry(S, i, i), fmpz_mat_entry(S, i + 1, i + 1)))
                continue;
            fmpz_xgcd(d, x, y, fmpz_mat_entry(S, i, i), fmpz_mat_entry(S, i + 1, i + 1));
            fmpz_divexact(fmpz_mat_entry(S, i + 1, i + 1), fmpz_mat_entry(S, i + 1, i + 1), d);

            /* record row operations */
            for (k = 0; k < A->r; k++)
            {
                fmpz_addmul(fmpz_mat_entry(U, i, k), y, fmpz_mat_entry(U, i + 1, k));
                fmpz_submul(fmpz_mat_entry(U, i + 1, k), fmpz_mat_entry(S, i + 1, i + 1),
                        fmpz_mat_entry(U, i, k));
                fmpz_neg(fmpz_mat_entry(U, i + 1, k), fmpz_mat_entry(U, i + 1, k));
            }

            fmpz_divexact(y, fmpz_mat_entry(S, i, i), d);

            /* record column operations */
            for (k = 0; k < A->c; k++)
            {
                fmpz_addmul(fmpz_mat_entry(V, k, i + 1), x, fmpz_mat_entry(V, k, i));
                fmpz_submul(fmpz_mat_entry(V, k, i), y, fmpz_mat_entry(V, k, i + 1));
                fmpz_swap(fmpz_mat_entry(V, k, i + 1), fmpz_mat_entry(V, k, i));
            }

            fmpz_mul(fmpz_mat_entry(S, i + 1, i + 1), fmpz_mat_entry(S, i + 1, i + 1), fmpz_mat_entry(S, i, i));
            fmpz_set(fmpz_mat_entry(S, i, i), d);
        }
    }
    fmpz_clear(d);
    fmpz_clear(x);
    fmpz_clear(y);
}

void fmpz_mat_snf_diagonal_trans(fmpz_mat_t S, fmpz_mat_t U, fmpz_mat_t V, const fmpz_mat_t A)
{
    fmpz_mat_one(U);
    fmpz_mat_one(V);
    _fmpz_mat_snf_diagonal_trans(S, U, V, A);
}
