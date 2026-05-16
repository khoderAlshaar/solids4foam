/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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
    ThormberRoe3ALEFlux

Description
    Implementation of inviscid numerical fluxes based on Roe & Pike FDS scheme
    with Harten's entropy fix. This is a first order scheme

Author
    Aleksandar Jemcov
    Rewrite by Hrvoje Jasak

SourceFiles
    ThormberRoe3ALEFlux.H

\*---------------------------------------------------------------------------*/

#ifndef ThormberRoe3ALEFlux_H
#define ThormberRoe3ALEFlux_H

#include "dbnsFlux.H" 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class ThormberRoe3ALEFlux Declaration
\*---------------------------------------------------------------------------*/

class ThormberRoe3ALEFlux:
public dbnsFlux
{
public:
    TypeName("ThormberRoe3ALEFlux");

    // Constructors

    //- Construct null
    ThormberRoe3ALEFlux(const fvMesh&, const dictionary&)
        {}
    
    
    // Destructor - default
    virtual ~ThormberRoe3ALEFlux() {};

    // Member Functions

        //- Update fluxes based on current left and right states
    virtual void evaluateFlux
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
        const scalar& RLeft,
        const scalar& RRight,
        const scalar& CvLeft,
        const scalar& CvRight,
        const vector& Sf,
        const scalar& magSf,
        const scalar& meshPhi,
        const tensor& R,
        const tensor& RTranspos,
        const scalar& fp1Left,
        const scalar& fp1Right
    ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
