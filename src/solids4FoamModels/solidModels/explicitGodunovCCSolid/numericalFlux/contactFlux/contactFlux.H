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
            const dynamicFvMesh& mesh_;
            
            //- Linear moemntum
            volVectorField& lm_;
            
            pointVectorField& lmN_;

            //- deformation tensor
            volTensorField& F_;

            //- 1st Piola-Kirchhoff stress tensor
            volTensorField& P_;

            solidMaterialModel& model_;
            operations& op_;
            mechanics& mech_;

            gradientSchemes& grad_;
            //-----------------------------
            // Material face area
            const surfaceScalarField& magSf_;

            // Material face area normal vector
            const surfaceVectorField& Sf_;

            // Material normals
            const surfaceVectorField N_;

            // Density
            const dimensionedScalar& rho_;
            
            // Gradient of cell linear momentum
            volTensorField lmGrad_;

            // First Piola Kirchhoff stress tensor
            volVectorField Px_;
            volVectorField Py_;
            volVectorField Pz_;

            // Gradients of first Piola Kirchhoff stress tensor
            volTensorField PxGrad_;
            volTensorField PyGrad_;
            volTensorField PzGrad_;

            //Creating fields for reconstruction
            // Reconstruction of linear momentum
            surfaceVectorField lm_M_;
            surfaceVectorField lm_P_;

            // Reconstruction of PK1 stresses
            surfaceTensorField P_M_;
            surfaceTensorField P_P_;

            // Reconstruction of traction
            surfaceVectorField t_M_;
            surfaceVectorField t_P_;

            // Creating fields for riemann solver
            // Surface tensors
            surfaceTensorField S_lm_;
            surfaceTensorField S_t_;

            // Contact traction
            surfaceVectorField tC_;

            // Contact linear momentum
            surfaceVectorField lmC_;

            // Volumetric fields
            volVectorField t_b_;
            volVectorField lm_b_;

            //limiter fields
            volVectorField phi_lm_;
            volTensorField phi_P_;

            // Creating fields for the constrained procedure
            // Constrained class
            interpolationSchemes interpolate_;

            // Cell-averaged linear momentum
            volVectorField lmR_;

            // Local gradient of cell-averaged linear momentum
            volTensorField lmRgrad_;

        // Fluxes

            //- Density flux
            surfaceVectorField lmFlux_;

            //- Velocity flux
            surfaceTensorField FFlux_;


    // Private Member Functions

       void  reconstruction();

       void curllFreeAlgorithm();

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
            const word& region,
            const dynamicFvMesh& mesh,
            volVectorField& lm,
            pointVectorField& lmN,
            volTensorField& F,
            volTensorField& P,
            solidMaterialModel& model,
            operations& op,
            mechanics& mech,
            gradientSchemes& grad
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
