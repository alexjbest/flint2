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

#ifndef NMOD_SPARSE_MAT_H
#define NMOD_SPARSE_MAT_H

#ifdef NMOD_SPARSE_MAT_INLINES_C
#define NMOD_SPARSE_MAT_INLINE FLINT_DLL
#else
#define NMOD_SPARSE_MAT_INLINE static __inline__
#endif

#undef ulong
#define ulong ulongxx /* interferes with system includes */
#include <stdlib.h>
#undef ulong
#include <gmp.h>
#define ulong mp_limb_t

#include "flint.h"
#include "longlong.h"
#include "ulong_extras.h"
#include "nmod_vec.h"

#ifdef __cplusplus
 extern "C" {
#endif

typedef struct
{
    slong pos;
    mp_limb_t entry;
}
nmod_sparse_mat_entry_struct;

typedef struct
{
    slong r;
    slong c;
    nmod_sparse_mat_entry_struct ** rows;
    mp_limb_t * row_supports;
    mp_limb_t * row_alloc;
    nmod_t mod;
}
nmod_sparse_mat_struct;

/* nmod_sparse_mat_t allows reference-like semantics for nmod_sparse_mat_struct */
typedef nmod_sparse_mat_struct nmod_sparse_mat_t[1];

NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_ensure_row_alloc(const nmod_sparse_mat_t mat, slong i, slong len)
{
    if (mat->row_alloc[i] >= len)
        return;
    else
    {
        mat->rows[i] = flint_realloc(mat->rows[i], (FLINT_MIN(len, mat->c)) * sizeof(nmod_sparse_mat_entry_struct));
        mat->row_alloc[i] = FLINT_MIN(len, mat->c);
    }
}

NMOD_SPARSE_MAT_INLINE
mp_limb_t nmod_sparse_mat_get_entry(const nmod_sparse_mat_t mat, slong i, slong j)
{
    slong k;
    for (k = 0; k < mat->row_alloc[i]; k++)
    {
        if (mat->rows[i][k].pos == j)
            return mat->rows[i][k].entry;
    }
    return 0;
}

NMOD_SPARSE_MAT_INLINE
void _nmod_sparse_mat_clear_entry(const nmod_sparse_mat_t mat, slong i, slong j, slong k)
{
    mat->row_supports[i]--;

    for (; k < mat->row_supports[i]; k++)
        mat->rows[i][k] = mat->rows[i][k + 1];
}

NMOD_SPARSE_MAT_INLINE
void _nmod_sparse_mat_set_entry(const nmod_sparse_mat_t mat, slong i, slong j, slong k, mp_limb_t n)
{
    nmod_sparse_mat_entry_struct prev;

    mat->row_supports[i]++;
    nmod_sparse_mat_ensure_row_alloc(mat, i, mat->row_supports[i]);

    prev.pos = j;
    prev.entry = n;

    for (; k < mat->row_supports[i]; k++)
    {
        const nmod_sparse_mat_entry_struct t = mat->rows[i][k];
        mat->rows[i][k] = prev;
        prev = t;
    }
}

NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_set_entry(const nmod_sparse_mat_t mat, slong i, slong j, mp_limb_t n)
{
    slong k;
    NMOD_RED(n, n, mat->mod);
    if (n == 0)
    {
        /* TODO */
    }

    for (k = 0; k < mat->row_supports[i]; k++)
    {
        if (mat->rows[i][k].pos == j)
        {
            mat->rows[i][k].entry = n;
            return;
        }
        else if (mat->rows[i][k].pos > j) /* insert new entry */
        {
            break;
        }
    }
    _nmod_sparse_mat_set_entry(mat, i, j, k, n);
}

NMOD_SPARSE_MAT_INLINE
slong nmod_sparse_mat_nrows(const nmod_sparse_mat_t mat)
{
   return mat->r;
}

NMOD_SPARSE_MAT_INLINE
slong nmod_sparse_mat_ncols(const nmod_sparse_mat_t mat)
{
   return mat->c;
}

NMOD_SPARSE_MAT_INLINE
void _nmod_sparse_mat_set_mod(nmod_sparse_mat_t mat, mp_limb_t n)
{
    mat->mod.n = n;
    count_leading_zeros(mat->mod.norm, n);
    invert_limb(mat->mod.ninv, n << mat->mod.norm);
}

/* Memory management */
FLINT_DLL void nmod_sparse_mat_init(nmod_sparse_mat_t mat, slong rows, slong cols, mp_limb_t n);
FLINT_DLL void nmod_sparse_mat_init_set(nmod_sparse_mat_t mat, const nmod_sparse_mat_t src);
FLINT_DLL void nmod_sparse_mat_clear(nmod_sparse_mat_t mat);
FLINT_DLL void nmod_sparse_mat_one(nmod_sparse_mat_t mat);
FLINT_DLL void nmod_sparse_mat_swap(nmod_sparse_mat_t mat1, nmod_sparse_mat_t mat2);

FLINT_DLL void nmod_sparse_mat_print_pretty(const nmod_sparse_mat_t mat);

FLINT_DLL int nmod_sparse_mat_equal(const nmod_sparse_mat_t mat1, const nmod_sparse_mat_t mat2);

FLINT_DLL void nmod_sparse_mat_zero(nmod_sparse_mat_t mat);

FLINT_DLL int nmod_sparse_mat_is_zero(const nmod_sparse_mat_t mat);

NMOD_SPARSE_MAT_INLINE
int
nmod_sparse_mat_is_zero_row(const nmod_sparse_mat_t mat, slong i)
{
    return mat->row_alloc[i] == 0;
}

NMOD_SPARSE_MAT_INLINE
int nmod_sparse_mat_is_empty(const nmod_sparse_mat_t mat)
{
    return (mat->r == 0) || (mat->c == 0);
}

NMOD_SPARSE_MAT_INLINE
int nmod_sparse_mat_is_square(const nmod_sparse_mat_t mat)
{
    return (mat->r == mat->c);
}


FLINT_DLL void nmod_sparse_mat_set(nmod_sparse_mat_t B, const nmod_sparse_mat_t A);
FLINT_DLL void nmod_sparse_mat_transpose(nmod_sparse_mat_t B, const nmod_sparse_mat_t A);

/* Addition and subtraction */

FLINT_DLL void nmod_sparse_mat_add(nmod_sparse_mat_t C, const nmod_sparse_mat_t A, const nmod_sparse_mat_t B);
FLINT_DLL void nmod_sparse_mat_sub(nmod_sparse_mat_t C, const nmod_sparse_mat_t A, const nmod_sparse_mat_t B);
FLINT_DLL void nmod_sparse_mat_neg(nmod_sparse_mat_t B, const nmod_sparse_mat_t A);

/* Matrix-scalar arithmetic */

FLINT_DLL void nmod_sparse_mat_scalar_mul(nmod_sparse_mat_t B, const nmod_sparse_mat_t A, mp_limb_t c);
FLINT_DLL void nmod_sparse_mat_scalar_mul_add(nmod_sparse_mat_t dest, const nmod_sparse_mat_t X,
                                const mp_limb_t b, const nmod_sparse_mat_t Y);

/* Matrix multiplication */

FLINT_DLL void nmod_sparse_mat_mul(nmod_sparse_mat_t C, const nmod_sparse_mat_t A, const nmod_sparse_mat_t B);

/* Permutations */

NMOD_SPARSE_MAT_INLINE
void nmod_sparse_mat_swap_rows(nmod_sparse_mat_t mat, slong * perm, slong r, slong s)
{
    if (r != s)
    {
        nmod_sparse_mat_entry_struct * u;
        mp_limb_t u2;
        slong t;

        if (perm)
        {
            t = perm[s];
            perm[s] = perm[r];
            perm[r] = t;
        }

        u = mat->rows[s];
        mat->rows[s] = mat->rows[r];
        mat->rows[r] = u;

        u2 = mat->row_alloc[s];
        mat->row_alloc[s] = mat->row_alloc[r];
        mat->row_alloc[r] = u2;

        u2 = mat->row_supports[s];
        mat->row_supports[s] = mat->row_supports[r];
        mat->row_supports[r] = u2;
    }
}

#ifdef __cplusplus
}
#endif

#endif

