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
    solidNumericFlux

Description
    Single level numeric flux class for density-based solvers

Author
    Aleksandar Jemcov
    Rewrite by Hrvoje Jasak

SourceFiles
    solidNumericFlux.H
    solidNumericFlux.C

\*---------------------------------------------------------------------------*/

#ifndef solidNumericFlux_H
#define solidNumericFlux_H

#include "solidNumericFluxBase.H"
#include "solidFlux.H"
#include "dynamicFvMesh.H"



#include "solidMaterialModel.H" //! to be changed with solidMaterialModel to work with solids4Foam
#include "mechanics.H"
#include "operations.H"
#include "gradientSchemes.H"
// #include "interpolationSchemes.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class solidNumericFlux Declaration
\*---------------------------------------------------------------------------*/

class solidNumericFlux
:
    public solidNumericFluxBase
{
    
    // Private data

            const dynamicFvMesh& mesh_;

            autoPtr<solidFlux> pFlux_;

            //- Reference to mesh

    
        
        // Reference to primitive fields

             volVectorField& lm_;
            // surfaceVectorField lm_M_;

             pointVectorField& lmN_;

             volTensorField& F_;

             volTensorField& P_;

            const solidMaterialModel& model_;

            const operations& op_;

             mechanics& mech_;
             gradientSchemes& grad_;

                // Material face area
            const surfaceScalarField& magSf_;

            // Material face area normal vector
            const surfaceVectorField& Sf_;

            // Material normals
            const surfaceVectorField N_;

            // Spatial normals
            surfaceVectorField n_;


            // surfaceVectorField t_M_;

            // // Density
            const dimensionedScalar& rho_;

            
            // Eigenvalues
            surfaceScalarField Up_;
            surfaceScalarField Us_;
            surfaceScalarField lambda_; // Lame's 1st parameter


            // Creating fields for riemann solver
            // Surface tensors
            surfaceTensorField S_lm_;
            surfaceTensorField S_t_;


            surfaceVectorField T1_;
        
            surfaceVectorField T2_;

                    
            surfaceTensorField R_;

            surfaceTensorField RTranspos_;

    // Face reconstructed fields (owned!)
        surfaceVectorField lm_P_;
        surfaceVectorField lm_M_;

        surfaceTensorField F_P_;
        surfaceTensorField F_M_;

        surfaceTensorField P_P_;
        surfaceTensorField P_M_;

        // Tractions
        surfaceVectorField t_P_;
        surfaceVectorField t_M_;

        // Boundary values
        volVectorField lm_b_;
        volVectorField t_b_;

        volVectorField phi_lm_;

        volTensorField phi_P_;
        volTensorField phi_F_;


        // Fluxes
            // - linear momentum flux
            surfaceVectorField lmFlux_;

            // - deformation tensor flux
            surfaceTensorField FFlux_;



    // Private Member Functions

        //- Disallow default bitwise copy construct
        solidNumericFlux(const solidNumericFlux&);


        //- Disallow default bitwise assignment
        void operator=(const solidNumericFlux&);




public:

    //- Runtime type information
    // TypeName("solidNumericFlux");

    // Constructors

        //- Construct from components
        solidNumericFlux
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


    //- Destructor
    virtual ~solidNumericFlux()
    {}


    // Member Functions
    // Private Member Functions
       void setRotationalMatrix();

        // Return fluxes

            // - Return linear momentum flux
            virtual const surfaceVectorField& lmFlux() const
            {
                return lmFlux_;
            }

            // - Return deformation tensor flux
            virtual const surfaceTensorField& FFlux() const
            {
                return FFlux_;
            }


        // Update fluxes based on current state

            //- Compute flux
            virtual void computeFlux();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// #ifdef NoRepository
// #   include "solidNumericFlux.C"
// #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
