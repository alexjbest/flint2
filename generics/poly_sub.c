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

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "generics.h"

void
_elem_poly_sub(elem_ptr res, elem_srcptr poly1, long len1,
    elem_srcptr poly2, long len2, const ring_t ring)
{
    long i, min;
    long size = ring->size;

    min = FLINT_MIN(len1, len2);

    for (i = 0; i < min; i++)
        elem_sub(INDEX(res, i, size), SRC_INDEX(poly1, i, size), SRC_INDEX(poly2, i, size), ring);

    if (poly1 != res)
        for (i = min; i < len1; i++)
            elem_set(INDEX(res, i, size), SRC_INDEX(poly1, i, size), ring);

    for (i = min; i < len2; i++)
        elem_neg(INDEX(res, i, size), SRC_INDEX(poly2, i, size), ring);
}
