/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     5.0
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

Application
    hllcSGFlux

Description
    HLLC flux class, computes numerical fluxes based on HLLC scheme. Uses
    primitive variables as input and gives back conservative numerical fluxes.
    
    Stiffened equation of state is used to compute the required variables from 
    given input.
    
    Toro, E.I. and Spruce, M. and Spears, W.: Restoration of the Contact
    Surface in the HLL-Riemann Solver, Shock Waves, 4:25-34, 1996
    
Author
    Ghanshyam Bharate, IIT Bombay

SourceFiles
    hllcSGFlux.H

\*---------------------------------------------------------------------------*/

#ifndef hllcSGFlux_H
#define hllcSGFlux_H

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class hllcFlux Declaration
\*---------------------------------------------------------------------------*/

class hllcSGFlux
{
public:

    // Constructors

        //- Constructor null
        hllcSGFlux()
        {}


    // Destructor - default


    // Member Functions

        //- Update fluxes based on current left and right states
       void evaluateFlux
         (
            scalar& rhoFlux,
            vector& rhoUFlux,
            scalar& rhoEFlux,
            const scalar& pLeft,
            const scalar& pRight,
            const vector& ULeft,
            const vector& URight,
            const scalar& TLeft,
            const scalar& TRight,
            const scalar& q,
            const scalar& pinf,
            const scalar& gamma,
            const scalar& Cv,
            const vector& Sf,
            const scalar& magSf
         ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
