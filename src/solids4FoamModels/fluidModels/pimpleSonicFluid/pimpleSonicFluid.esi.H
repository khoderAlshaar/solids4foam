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
    pimpleSonicFluid

Description
    pimpleFoam fluid model from OpenFOAM-v2012.

Author
    Ported by Philip Cardiff, UCD.
    Extension to Robin coupling and thermal analysed by Zeljko Tukovic, Zagreb.

SourceFiles
    pimpleSonicFluid.esi.C

\*---------------------------------------------------------------------------*/

#ifndef pimpleSonicFluid_H
#define pimpleSonicFluid_H

#include "fluidModel.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "CorrectPhi.H"

#include "psiThermo.H"
#include "turbulentFluidThermoModel.H"
#include "CorrectPhi.H"
#include "constrainHbyA.H"
#include "fvCFD.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidModels
{

/*---------------------------------------------------------------------------*\
                        Class pimpleSonicFluid Declaration
\*---------------------------------------------------------------------------*/

class pimpleSonicFluid
:
    public fluidModel
{
    // Private data

        autoPtr<psiThermo> pThermo_;

        psiThermo& thermo_;
        
        //- Turbulence model
        autoPtr<compressible::turbulenceModel> turbulence_;

        volScalarField K_;

        volScalarField& e_;

        const volScalarField& psi_;

        //- Dynamic viscosity
        volScalarField mu_;

        //- Fluid density
        volScalarField rho_;


        //- Pointer to the velocity surface field
        autoPtr<surfaceVectorField> rhoUf_;

        //- rAU field
        volScalarField rAU_;

        //- Correct phi flag
        const bool correctPhi_;

        //- Check Courant number flag
        const bool checkMeshCourantNo_;

        // Move mesh in the outer corrector loop flag
        const bool moveMeshOuterCorrectors_;

        //- Cumulative continuity errors
        scalar cumulativeContErr_;

                //- Adjustable time-step switch
        const Switch adjustTimeStep_;

                //- For adjustable time-step, this is the maximum Courant number
        const scalar maxCo_;

        //- For adjustable time-step, this is the maximum time-stpe
        const scalar maxDeltaT_;
        



    // Private Member Functions
        //- Calculates and prints the continuity errors
        void compressibleContinuityErrs();

          //- Phi correction for mesh motion
        void CorrectFlux();

        //- Print compressible Courant number
        void compressibleCourantNo();

        //- Disallow default bitwise copy construct
        pimpleSonicFluid(const pimpleSonicFluid&);

        //- Disallow default bitwise assignment
        void operator=(const pimpleSonicFluid&);


public:

    //- Runtime type information
    TypeName("pimpleSonicFluid");

    // Constructors

        //- Construct from components
        pimpleSonicFluid
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Destructor

        virtual ~pimpleSonicFluid()
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
