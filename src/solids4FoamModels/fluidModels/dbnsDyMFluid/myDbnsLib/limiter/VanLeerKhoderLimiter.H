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
    VanLeerKhoderLimiter

Description
    van Leer limiter for multidimensional slope limiting
    Formula: phi(r) = 2*(r + |r|) / (1 + |r|)^2

\*---------------------------------------------------------------------------*/

#ifndef VanLeerKhoderLimiter_H
#define VanLeerKhoderLimiter_H

#include "scalar.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class VanLeerKhoderLimiter Declaration
\*---------------------------------------------------------------------------*/

class VanLeerKhoderLimiter
{
public:

    // Constructor
    VanLeerKhoderLimiter()
    {}

    // Limiter function for scalar fields
    inline void limiter
    (
        scalar& lim,
        const scalar& cellVolume,
        const scalar& deltaOneMax,
        const scalar& deltaOneMin,
        const scalar& extrapolate
    )
    {
        scalar r = 1;

            if (mag(extrapolate) < SMALL)
            {
                return;
            }

            if (extrapolate > SMALL)
            {
                r =  deltaOneMax/extrapolate;
            }
            else if (extrapolate < VSMALL)
            {
                r =  deltaOneMin/extrapolate;
            }
 
        scalar phi_r = (r + mag(r))/(1 + mag(r));

        lim = min(lim, phi_r);
    }

    // Limiter function for vector and tensor fields
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
