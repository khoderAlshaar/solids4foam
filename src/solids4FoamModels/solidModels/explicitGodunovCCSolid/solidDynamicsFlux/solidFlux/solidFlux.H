/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: Open Source CFD
   \\    /   O peration     | 
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  | 
-------------------------------------------------------------------------------
License
    This file isn't part of foam-extend nor OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    solidFlux

Description
    Basic class for of inviscid numerical fluxes.

Author
    Jiri Furst

SourceFiles
    solidFlux.H

\*---------------------------------------------------------------------------*/

#ifndef solidFlux_H
#define solidFlux_H

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class roeFlux Declaration
\*---------------------------------------------------------------------------*/

class solidFlux
{
private:
    solidFlux(const solidFlux&);
    void operator=(const solidFlux&);

public:
    TypeName("solidFlux");

    declareRunTimeSelectionTable
    (
        autoPtr,
        solidFlux,
        dictionary,
        (const fvMesh& mesh, const dictionary& dict),
        (mesh, dict)
    );

    // Constructors
    solidFlux() {};

    //- Return a reference to the selected solidFlux
    static autoPtr<solidFlux> New(
        const fvMesh& mesh, 
        const dictionary& dict,
	const word key="solidFlux"
    );

    // Destructor - default
    virtual ~solidFlux() {};

    // Member Functions

        // Update fluxes based on current left and right states
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
        ) const = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
