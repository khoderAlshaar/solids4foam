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
    sonicFluid

Description
    Transient solver for trans-sonic/supersonic, laminar flow of a
    compressible liquid. Assumes a linear equation of state between
    density and pressure related by the bulk modulus.

Author
    Iago Lessa de Oliveira, FEIS/UNESP/UCD.
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    sonicFluid.esi.C

\*---------------------------------------------------------------------------*/

#ifndef sonicFluid_H
#define sonicFluid_H

#include "fluidModel.H"
#include "pimpleControl.H"
#include "uniformDimensionedFields.H"
#include "pressureControl.H"
#include "fvcSmooth.H"

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
                        Class sonicFluid Declaration
\*---------------------------------------------------------------------------*/

class sonicFluid
:
    public fluidModel
{
    // Private data
        autoPtr<psiThermo> pThermo_;

        psiThermo& thermo_;
        
        autoPtr<compressible::turbulenceModel> turbulence_;

        volScalarField K_;

        volScalarField& e_;

        const volScalarField& psi_;

        // Transport properties dictionary
        // IOdictionary transportProperties_;

        //- Dynamic viscosity
        volScalarField mu_;




        //- Fluid density
        volScalarField rho_;

        //- Creates and initialises the momentum field rhoUf
        autoPtr<surfaceVectorField> rhoUf_;

        IOMRFZoneList MRF;


        //- rAU field
        // volScalarField rAU_;

        //- Cumulative continuity error
        scalar cumulativeContErr_;

    // Private Member Functions

        //- Calculates and prints the continuity errors
        void compressibleContinuityErrs();

        //- Solve density equation
        void solveRhoEqn();

        //- Phi correction for mesh motion
        void CorrectFlux();

        //- Print compressible Courant number
        void compressibleCourantNo();

        //- Solve pressure equation
        void solvePEqn
        (
            const fvVectorMatrix& UEqn
        );

        //- Disallow default bitwise copy construct
        sonicFluid(const sonicFluid&);

        //- Disallow default bitwise assignment
        void operator=(const sonicFluid&);


public:

    //- Runtime type information
    TypeName("sonicFluid");

    // Constructors

        //- Construct from components
        sonicFluid
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Destructor

        virtual ~sonicFluid()
        {}


    // Member Functions

        // Access

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
