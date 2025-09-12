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
    dbnsDyMFluid

Description
    pimpleFoam fluid model from OpenFOAM-v2012.

Author
    Ported by Philip Cardiff, UCD.
    Extension to Robin coupling and thermal analysed by Zeljko Tukovic, Zagreb.

SourceFiles
    dbnsDyMFluid.esi.C

\*---------------------------------------------------------------------------*/

#ifndef dbnsDyMFluid_H
#define dbnsDyMFluid_H

#include "fluidModel.H"
#include "fvCFD.H"
#include "basicPsiThermo.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"
#include "myTurbulenceModel.H"
// #include "bound.H"
// #include "roeALEFlux.H"
// #include "MDLimiter.H"
// #include "firstOrderLimiter.H"
// #include "BarthJespersenLimiter.H"
// #include "VenkatakrishnanLimiter.H"
#include "numericFlux.H"
#include "backwardDualDdtScheme.H"
#include "myBackwardDualDdtScheme.H"  
// #include "dynamicFvMesh.H"

#include "localTimeStep.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidModels
{

/*---------------------------------------------------------------------------*\
                        Class dbnsDyMFluid Declaration
\*---------------------------------------------------------------------------*/

class dbnsDyMFluid
:
    public fluidModel
{
    // Private data
        Time& runTime_;
        autoPtr<basicPsiThermo> pThermo_;

        basicPsiThermo& thermo_;
        
        volScalarField& h_;
        // volScalarField& p_;
        const volScalarField& T_;

                // //- Fluid density
        volScalarField rho_;
        // volVectorField U_;
        volVectorField rhoU_;
        volScalarField rhoE_;



        // autoPtr<basicNumericFlux> dbnsFluxPtr_;
    
        // basicNumericFlux& dbnsFlux_;

        autoPtr<numericFlux> dbnsFluxPtr_;
    
        numericFlux& dbnsFlux_;



        // //- Turbulence model
        autoPtr<compressible::myTurbulenceModel> turbulence_;

            IOField<scalar> physDeltaT_;
        
        localTimeStep localTimeStep_;
        volScalarField CoDeltaT_;
        dimensionedScalar pseudoTimeStep_;
        label numberSubCycles_;
        scalar tolerance_;
        scalar relTol_;


        // volScalarField K_;

        // volScalarField& e_;

        // const volScalarField& psi_;

        // //- Dynamic viscosity
        // volScalarField mu_;


        // //- Pointer to the velocity surface field
        // autoPtr<surfaceVectorField> rhoUf_;


                //- Adjustable time-step switch
         Switch adjustTimeStep_;

                //- For adjustable time-step, this is the maximum Courant number
         scalar maxCo_;

        //- For adjustable time-step, this is the maximum time-stpe
        const scalar maxDeltaT_;
        



    // Private Member Functions
    bool converged
(
    const int iCorr,
    const dimensionedScalar pDeltaT,
    // const volVectorField& vf
    const volScalarField& vf,
    const label numberSubCycles,
    scalar tolerance   
);
        // //- Calculates and prints the continuity errors
        // void compressibleContinuityErrs();

        //   //- Phi correction for mesh motion
        // void CorrectFlux();

        // //- Print compressible Courant number
        // void compressibleCourantNo();

        //- Disallow default bitwise copy construct
        dbnsDyMFluid(const dbnsDyMFluid&);

        //- Disallow default bitwise assignment
        void operator=(const dbnsDyMFluid&);


public:

    //- Runtime type information
    TypeName("dbnsDyMFluid");

    // Constructors

        //- Construct from components
        dbnsDyMFluid
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Destructor

        virtual ~dbnsDyMFluid()
        {}


    // Member Functions

        // Access
            virtual void setDeltaT(Time& runTime);

            //- Patch viscous force (N/m2)
            virtual tmp<vectorField> patchViscousForce
            (
                const label patchID
            ) const;

            //- Patch pressure force (N/m2)
            virtual tmp<scalarField> patchPressureForce
            (
                const label patchID
            ) const;

 
        // Edit

            //- Evolve the fluid model
            virtual bool evolve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
