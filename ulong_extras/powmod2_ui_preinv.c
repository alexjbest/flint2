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

    Copyright (C) 2009, 2013, 2016 William Hart

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

ulong
n_powmod2_ui_preinv(ulong a, ulong exp, ulong n, ulong ninv)
{
    ulong x, norm;

    FLINT_ASSERT(n != 0);

    if (exp == 0)
    {
        /* anything modulo 1 is 0 */
        return n == 1 ? 0 : 1;
    }

    if (a == 0)
        return 0;

    if (a >= n)
        a = n_mod2_preinv(a, n, ninv);

    count_leading_zeros(norm, n);
    a <<= norm;
    n <<= norm;

    while ((exp & 1) == 0)
    {
        a = n_mulmod_preinv(a, a, n, ninv, norm);
        exp >>= 1;
    }

    x = a;

    while (exp >>= 1)
    {
        a = n_mulmod_preinv(a, a, n, ninv, norm);

        if (exp & 1)
            x = n_mulmod_preinv(x, a, n, ninv, norm);
    }

    return x >> norm;
}
