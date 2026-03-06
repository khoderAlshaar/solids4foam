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
    VanLeerLimiter

Description
    van Leer limiter for multidimensional slope limiting
    Formula: phi(r) = 2*(r + |r|) / (1 + |r|)^2

\*---------------------------------------------------------------------------*/

#ifndef VanLeerLimiter_H
#define VanLeerLimiter_H

#include "scalar.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class VanLeerLimiter Declaration
\*---------------------------------------------------------------------------*/

class VanLeerLimiter
{
public:

    // Constructor
    VanLeerLimiter()
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
        // Compute gradient ratio
        scalar deltaMean = 0.5*(deltaOneMax - deltaOneMin);
        
        if (mag(deltaMean) < SMALL)
        {
            lim = 0.0;
            return;
        }
        
        scalar r = extrapolate / (deltaMean + SMALL);
        
        // van Leer formula: 2*(r + |r|) / (1 + |r|)^2
        scalar phi_r = (2.0*(r + mag(r))) / (sqr(1.0 + mag(r)) + SMALL);
        
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
