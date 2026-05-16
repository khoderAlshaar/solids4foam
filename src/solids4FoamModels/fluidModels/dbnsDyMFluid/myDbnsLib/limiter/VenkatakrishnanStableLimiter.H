/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    VenkatakrishnanStableLimiter

Description
    Venkatakrishnan differentiable limiter

Author
    Aleksandar Jemcov
    Rewrite by Hrvoje Jasak

\*---------------------------------------------------------------------------*/

#ifndef VenkatakrishnanStableLimiter_H
#define VenkatakrishnanStableLimiter_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class VenkatakrishnanStableLimiter Declaration
\*---------------------------------------------------------------------------*/

class VenkatakrishnanStableLimiter
{
    // Private data

        //- Limiter value
        // scalar k_;


public:

    // Constructor

        //- Construct null
        VenkatakrishnanStableLimiter()
        // :
        //     k_(5)
            // k_(0.001)
        {}

        inline void limiter
        (
            scalar& lim,
            const scalar& cellVolume,
            const scalar& deltaOneMax,
            const scalar& deltaOneMin,
            const scalar& extrapolate
        )
        {
            if (mag(extrapolate) < SMALL)
            {
                return;
            }

            else if (extrapolate > VSMALL)
            {
                scalar r =  deltaOneMax/extrapolate;

                scalar y = (sqr(r)+ 2*r)/ stabilise(sqr(r) + r +2,SMALL);

                lim = min(lim,y);
            }
            else if (extrapolate < VSMALL)
            {
                scalar r =  deltaOneMin/extrapolate;

                scalar y = (sqr(r)+ 2*r)/ stabilise(sqr(r) + r +2,SMALL);

                lim = min(lim, y);
            }
        }

        //- Set Type limiter
        template<class Type>
        inline void limiter
        (
            Type& lim,
            const scalar& cellVolume,
            const Type& deltaOneMax,
            const Type& deltaOneMin,
            const Type& extrapolate
        )
        {
            for (direction cmpt = 0; cmpt < Type::nComponents; cmpt++)
            {
                limiter
                (
                    lim.component(cmpt),
                    cellVolume,
                    deltaOneMax.component(cmpt),
                    deltaOneMin.component(cmpt),
                    extrapolate.component(cmpt)
                );
            }
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
