/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mod_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("resultant....");
    fflush(stdout);

    /* Check res(f, g) == (-1)^(deg f deg g) res(g, f) */
    for (i = 0; i < 60 * flint_test_multiplier(); i++)
    {
        fmpz_t x, y, n;
        fmpz_mod_poly_t f, g;

        fmpz_init(n);
        fmpz_init(x);
        fmpz_init(y);

        fmpz_set_ui(n, n_randtest_prime(state, 0));

        fmpz_mod_poly_init(f, n);
        fmpz_mod_poly_init(g, n);
        
        fmpz_mod_poly_randtest(f, state, n_randint(state, 300));
        fmpz_mod_poly_randtest(g, state, n_randint(state, 300));

        fmpz_mod_poly_resultant(x, f, g);
        fmpz_mod_poly_resultant(y, g, f);

        if ((fmpz_mod_poly_degree(f) * fmpz_mod_poly_degree(g)) % 2)
           fmpz_negmod(y, y, &f->p);

        result = (fmpz_equal(x, y));
        if (!result)
        {
            flint_printf("FAIL (res(f, g) == (-1)^(deg f deg g) res(g, f)):\n");
            fmpz_mod_poly_print(f), flint_printf("\n\n");
            fmpz_mod_poly_print(g), flint_printf("\n\n");
            printf("x = "); fmpz_print(x); printf("\n");
            printf("y = "); fmpz_print(y); printf("\n");
            printf("n = "); fmpz_print(n); printf("\n");
            abort();
        }

        fmpz_mod_poly_clear(f);
        fmpz_mod_poly_clear(g);
        fmpz_clear(n);
        fmpz_clear(x);
        fmpz_clear(y);
    }

    /* Check res(f h, g) == res(f, g) res(h, g) */
    for (i = 0; i < 60 * flint_test_multiplier(); i++)
    {
        fmpz_t x, y, z, n;
        fmpz_mod_poly_t f, g, h;

        fmpz_init(n);
        fmpz_init(x);
        fmpz_init(y);
        fmpz_init(z);

        fmpz_set_ui(n, n_randtest_prime(state, 0));

        fmpz_mod_poly_init(f, n);
        fmpz_mod_poly_init(g, n);
        fmpz_mod_poly_init(h, n);
        
        fmpz_mod_poly_randtest(f, state, n_randint(state, 200));
        fmpz_mod_poly_randtest(g, state, n_randint(state, 200));
        fmpz_mod_poly_randtest(h, state, n_randint(state, 100));

        fmpz_mod_poly_resultant(y, f, g);
        fmpz_mod_poly_resultant(z, h, g);
        fmpz_mul(y, y, z);
        fmpz_mod(y, y, n);
        fmpz_mod_poly_mul(f, f, h);
        fmpz_mod_poly_resultant(x, f, g);
        
        result = (fmpz_equal(x, y));
        if (!result)
        {
            flint_printf("FAIL (res(f h, g) == res(f, g) res(h, g)):\n");
            fmpz_mod_poly_print(f), flint_printf("\n\n");
            fmpz_mod_poly_print(g), flint_printf("\n\n");
            fmpz_mod_poly_print(h), flint_printf("\n\n");
            printf("x = "); fmpz_print(x); printf("\n");
            printf("y = "); fmpz_print(y); printf("\n");
            printf("z = "); fmpz_print(z); printf("\n");
            printf("n = "); fmpz_print(n); printf("\n");
            abort();
        }

        fmpz_mod_poly_clear(f);
        fmpz_mod_poly_clear(g);
        fmpz_mod_poly_clear(h);
        fmpz_clear(n);
        fmpz_clear(x);
        fmpz_clear(y);
        fmpz_clear(z);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

