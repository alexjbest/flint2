/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"


void mpoly_monomials_inflate(ulong * Aexps, mp_bitcnt_t Abits,
                       const ulong * Bexps, mp_bitcnt_t Bbits, slong Blength,
               const fmpz * shift, const fmpz * stride, const mpoly_ctx_t mctx)
{
    slong i, j;
    slong NA, NB;
    slong nvars = mctx->nvars;
    fmpz * exps;
    TMP_INIT;

    TMP_START;
    exps = (fmpz *) TMP_ALLOC(nvars*sizeof(fmpz));
    for (j = 0; j < nvars; j++)
        fmpz_init(exps + j);

    NA = mpoly_words_per_exp(Abits, mctx);
    NB = mpoly_words_per_exp(Bbits, mctx);

    for (i = 0; i < Blength; i++)
    {
        mpoly_get_monomial_ffmpz(exps, Bexps + NB*i, Bbits, mctx);
        for (j = 0; j < nvars; j++)
        {
            fmpz_mul(exps + j, exps + j, stride + j);
            fmpz_add(exps + j, exps + j, shift + j);
        }
        FLINT_ASSERT(Abits >= mpoly_exp_bits_required_ffmpz(exps, mctx));

        mpoly_set_monomial_ffmpz(Aexps + NA*i, exps, Abits, mctx);
    }

    for (j = 0; j < nvars; j++)
        fmpz_clear(exps + j);

    TMP_END;
}

