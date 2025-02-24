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
    roeFlux

Description
    Virtual base class for solid mechanics models

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    roeFlux.C
    roeFluxTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef roeFlux_H
#define roeFlux_H

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
                           Class roeFlux Declaration
\*---------------------------------------------------------------------------*/

class roeFlux
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


        // Roe specific fields
        //decompoz of tensor F
        volVectorField Fx_;
        volVectorField Fy_;
        volVectorField Fz_;

        // compute F gradient
        volTensorField FxGrad_;
        volTensorField FyGrad_;
        volTensorField FzGrad_;

        volTensorField phi_F_;

        // Reconstruction of F tensor
        surfaceTensorField F_M_;
        surfaceTensorField F_P_;

       

        surfaceVectorField T1_;
        
        surfaceVectorField T2_;

                
        surfaceTensorField R_;

        surfaceTensorField RTranspos_;

        //------- Surface Vector Fields ---------
        surfaceVectorField lm_M_hat_;
        surfaceVectorField lm_P_hat_;

        //------- Surface Tensor Fields ---------
        surfaceTensorField F_M_hat_;
        surfaceTensorField F_P_hat_;

        // Reconstruction of PK1 stresses
        surfaceTensorField P_hat_P_;
        surfaceTensorField P_hat_M_;

        // Define the fluxes after rotation on the left and right
        surfaceVectorField flux_lm_M_hat_;
        surfaceVectorField flux_lm_P_hat_;
        surfaceTensorField flux_F_M_hat_;
        surfaceTensorField flux_F_P_hat_;

        //------- Roe average values -----------
        // Difference fields
        surfaceVectorField deltaF_lm_;
        surfaceTensorField deltaF_F_;

        // Eigenvalues
        surfaceScalarField Up_;
        surfaceScalarField Us_;
        surfaceScalarField lambda_; // Lame's 1st parameter

        // Eigenvector linear momentum
        surfaceVectorField l1_lm_;
        surfaceVectorField l2_lm_;
        surfaceVectorField l3_lm_;
        surfaceVectorField l4_lm_;
        surfaceVectorField l5_lm_;
        surfaceVectorField l6_lm_;

        // Eigenvector deformation tensor presented in tensor form
        surfaceTensorField l1_F_;
        surfaceTensorField l2_F_;
        surfaceTensorField l3_F_;
        surfaceTensorField l4_F_;
        surfaceTensorField l5_F_;
        surfaceTensorField l6_F_;

        // Wave strength
        surfaceScalarField alpha1_;
        surfaceScalarField alpha2_;
        surfaceScalarField alpha3_;
        surfaceScalarField alpha4_;
        surfaceScalarField alpha5_;
        surfaceScalarField alpha6_;

        // Sum
        surfaceVectorField sum_lm_;
        surfaceTensorField sum_F_;

        // Roe Fluxes
        surfaceVectorField roeFlux_lm_;
        surfaceTensorField roeFlux_F_;
        surfaceVectorField roeFlux_lm_hat_;
        surfaceTensorField roeFlux_F_hat_;

    // Private Member Functions
       void setRotationalMatrix();

       void  reconstruction();

       void curllFreeAlgorithm();

        //- Disallow default bitwise copy construct
        roeFlux(const roeFlux&);

        //- Disallow default bitwise assignment
        void operator=(const roeFlux&);


protected:

    // Protected member functions


public:

    //- Runtime type information
    TypeName("roe");


    // Constructors

        //- Construct from components
        roeFlux
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

        virtual ~roeFlux();


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
