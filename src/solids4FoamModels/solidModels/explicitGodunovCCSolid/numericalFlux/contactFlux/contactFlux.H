/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

Class
    contactFlux

Description
    Virtual base class for solid mechanics models

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    contactFlux.C
    contactFluxTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef contactFlux_H
#define contactFlux_H

#include "numericalFlux.H"
// #include "IOdictionary.H"
// #include "scalarIOField.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{

namespace numericalFluxs
{
/*---------------------------------------------------------------------------* \
                           Class contactFlux Declaration
\*---------------------------------------------------------------------------*/

class contactFlux
:
    public numericalFlux
{
    // Private data
            //- Linear moemntum
            const volVectorField& lm_;

            //- deformation tensor
            const volTensorField& F_;

            //- 1st Piola-Kirchhoff stress tensor
            const volTensorField& P_;

            solidMaterialModel& model_;
            operations& op_;
            mechanics& mech_;

        // Fluxes

            //- Density flux
            surfaceVectorField lmFlux_;

            //- Velocity flux
            surfaceTensorField FFlux_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        contactFlux(const contactFlux&);

        //- Disallow default bitwise assignment
        void operator=(const contactFlux&);


protected:

    // Protected member functions


public:

    //- Runtime type information
    TypeName("contact");


    // Constructors

        //- Construct from components
        contactFlux
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

    // Destructor

        virtual ~contactFlux();


    // Member Functions

        // Return fluxes

            //- Return density flux
            virtual const surfaceVectorField& lmFlux() const
            {
                return lmFlux_;
            }

            //- Return velocity flux
            virtual const surfaceTensorField& FFlux() const
            {
                return FFlux_;
            }


       // Return residuals

            //- Return density equation residual
            virtual tmp<vectorField> lmResidual() const
            {
                return fvc::div(lmFlux_)().internalField();
            }

            //- Return momentum equation flux
            virtual tmp<tensorField> FResidual() const
            {
                return fvc::div(FFlux_)().internalField();
            }



        // Update fluxes based on current state

            //- Compute flux
            virtual void computeFlux();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace numericalFluxs


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
