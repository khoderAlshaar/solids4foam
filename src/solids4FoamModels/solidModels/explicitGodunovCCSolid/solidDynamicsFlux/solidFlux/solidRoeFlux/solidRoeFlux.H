/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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
    solidRoeFlux

Description
    Implementation of inviscid numerical fluxes based on Roe & Pike FDS scheme
    with Harten's entropy fix. This is a first order scheme

Author
    Aleksandar Jemcov
    Rewrite by Hrvoje Jasak

SourceFiles
    solidRoeFlux.H

\*---------------------------------------------------------------------------*/

#ifndef solidRoeFlux_H
#define solidRoeFlux_H

// #include "fvCFD.H"
#include "solidFlux.H" 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class solidRoeFlux Declaration
\*---------------------------------------------------------------------------*/

class solidRoeFlux:
public solidFlux
{
public:
    TypeName("solidRoeFlux");
    // Constructors

    //- Construct null
    solidRoeFlux(const fvMesh&, const dictionary&)
        {}

    // Destructor - default
    virtual ~solidRoeFlux() {};

    // Member Functions

        //- Update fluxes based on current left and right states
    virtual void evaluateFlux
        (
            vector& lmFlux,
            tensor& FFlux,
            const vector& lm_P,
            const vector& lm_M,
            const tensor& F_P,
            const tensor& F_M,
            const tensor& P_P,
            const tensor& P_M,
            // const vector& t_P,
            // const vector& t_M,
            const tensor& R,
            const tensor& RTranspos,
            // const tensor& S_lm,
            // const tensor& S_t,
            const scalar& rho,
            const scalar& lambda,
            const scalar& Up,
            const scalar& Us//,
            // const scalar& magSf,
            // const vector& Sf,
            // const vector& N
        ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
