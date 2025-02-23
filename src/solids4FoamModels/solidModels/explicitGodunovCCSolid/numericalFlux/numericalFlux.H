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
    Foam::numericalFlux

Description
    Base class for runtime-selectable flux and limiters

Author
    Henrik Rusche

SourceFiles
    numericalFlux.C
    newnumericalFlux.C

\*---------------------------------------------------------------------------*/

#ifndef numericalFlux_H
#define numericalFlux_H

#include "fvMesh.H"
#include "dynamicFvMesh.H"
#include "volFields.H"
#include "runTimeSelectionTables.H"
#include "basicThermo.H"


#include "solidMaterialModel.H" //! to be changed with solidMaterialModel to work with solids4Foam
#include "mechanics.H"
#include "operations.H"
#include "gradientSchemes.H"
#include "interpolationSchemes.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class numericalFlux Declaration
\*---------------------------------------------------------------------------*/

class numericalFlux
{
    // Private data

        //- Reference to mesh
        // const fvMesh& mesh_;


protected:

    // Protected member functions

        //- Construct as copy (not implemented)
        numericalFlux(const numericalFlux&);


public:

    //- Runtime type information
    TypeName("numericalFlux");


    //- Declare run-time constructor selection table

#ifndef SWIG
    declareRunTimeSelectionTable
    (
        autoPtr,
        numericalFlux,
        state,
        (
            Time& runTime,
            const dynamicFvMesh& mesh,
            const volVectorField& lm,
            const volTensorField& F,
            const volTensorField& P,
            solidMaterialModel& model,
            operations& op,
            mechanics& mech
        ),
        (runTime, mesh, lm, F, P, model, op, mech)
    );
#endif


    // Constructors

        //- Construct from components
        numericalFlux
        (
            Time& runTime
        );


    //- Selector
    static autoPtr<numericalFlux> New
    (
        Time& runTime,
        const dynamicFvMesh& mesh,
        const volVectorField& lm,
        const volTensorField& F,
        const volTensorField& P,
        solidMaterialModel& model,
        operations& op,
        mechanics& mech
    );

    //- Destructor
    virtual ~numericalFlux();



        // Return fluxes

            //- Return linear momentum  flux
            virtual const surfaceVectorField& lmFlux() const = 0;

            //- Return deformation tensor flux
            virtual const surfaceTensorField& FFlux() const = 0;


       // Return residuals

            //- Return linear momentum equation residual
            virtual tmp<vectorField> lmResidual() const = 0;

            //- Return deformation tensor equation flux
            virtual tmp<tensorField> FResidual() const = 0;


        // Update fluxes based on current state

            //- Compute flux
            virtual void computeFlux() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
