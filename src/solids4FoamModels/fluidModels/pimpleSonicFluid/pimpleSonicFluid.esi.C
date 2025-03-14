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

\*---------------------------------------------------------------------------*/

#include "pimpleSonicFluid.H"
#include "addToRunTimeSelectionTable.H"
#include "CorrectPhi.H"
#include "fvc.H"
#include "fvm.H"
#include "constrainHbyA.H"
#include "constrainPressure.H"
#include "findRefCell.H"
#include "elasticSlipWallVelocityFvPatchVectorField.H"
#include "elasticWallVelocityFvPatchVectorField.H"
#include "elasticWallPressureFvPatchScalarField.H"
#include "movingWallPressureFvPatchScalarField.H"
#include "EulerDdtScheme.H"
#include "backwardDdtScheme.H"
#include "thermalRobinFvPatchScalarField.H"

#include "adjustPhi.H"
#include "fvc.H"
#include "fvm.H"
#include "CorrectPhi.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pimpleSonicFluid, 0);
addToRunTimeSelectionTable(fluidModel, pimpleSonicFluid, dictionary);

// * * * * * * * * * * * * * * * Private Members * * * * * * * * * * * * * * //

void pimpleSonicFluid::CorrectFlux()
{
    // Store divrhoU from the previous mesh so that it can be mapped
    // and used in correctPhi to ensure the corrected phi has the
    // same divergence
    autoPtr<volScalarField> divrhoU;

    divrhoU.set
    (
        new volScalarField
        (
                "divrhoU",
                fvc::div(fvc::absolute(phi(), rho_, U()))
        )
    );

    // Store momentum to set rhoUf for introduced faces.
    autoPtr<volVectorField> rhoU;

    // Calculate absolute flux
    // from the mapped surface velocity
    phi() = mesh().Sf() & rhoUf_();

    // Define volScalarField to hold phi
    // to pass to compressible CorrectPhi function


    CorrectPhi
    (
        U(),
        phi(),
        p(),
        rho_,
        psi_,
        dimensionedScalar("rAUf", dimTime, 1),
        divrhoU(),
        pimple()
        //,
        //true
    );
}


// void pimpleSonicFluid::compressibleContinuityErrs()
// {
//     scalar sumLocalContErr =
//                 (
//                     sum
//                     (
//                         mag(rho_ - rho0_ - psi_*(p() - p0_))
//                     )/sum(rho_)
//                 ).value();

//     scalar globalContErr =
//                 (
//                         sum(rho_ - rho0_ - psi_*(p() - p0_))/sum(rho_)
//                 ).value();

//     cumulativeContErr_ += globalContErr;

//     Info<< "time step continuity errors : sum local = "
//                 << sumLocalContErr
//         << ", global = " << globalContErr
//         << ", cumulative = " << cumulativeContErr_ << endl;
// }



void pimpleSonicFluid::compressibleCourantNo()
{
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;

    surfaceScalarField SfUfbyDelta
        (
            mesh().surfaceInterpolation::deltaCoeffs()*mag(phi())
            /fvc::interpolate(rho_)
        );

    CoNum = max
            (
                SfUfbyDelta/mesh().magSf()
            ).value()*runTime().deltaT().value();

    meanCoNum = (
                    sum(SfUfbyDelta)/sum(mesh().magSf())
                ).value()*runTime().deltaT().value();

    Info<< "Region: " << mesh().name()
        << " Courant Number mean: " << meanCoNum
        << " max: " << CoNum << endl;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pimpleSonicFluid::pimpleSonicFluid
(
    Time& runTime,
    const word& region
)
:
    fluidModel(typeName, runTime, region),
    pThermo_
    (
        psiThermo::New(mesh())
    ),
     thermo_(pThermo_()),
    turbulence_
    (
        compressible::turbulenceModel::New
        (
            rho_,
            U(),
            phi(),
            thermo_
        )
    ),
    K_
    (
        "K",
        0.5*magSqr(U())
    ),

    e_(thermo_.he()),

    psi_(thermo_.psi()),

    
    mu_(thermo_.mu()),


    rho_
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        thermo_.rho()
    ),

    rhoUf_(),
    rAU_
    (
        IOobject
        (
            "rAU",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        runTime.deltaT(),
        calculatedFvPatchScalarField::typeName
    ),

    correctPhi_(pimple().dict().lookupOrDefault("correctPhi", false)),
    checkMeshCourantNo_
    (
        pimple().dict().lookupOrDefault("checkMeshCourantNo", false)
    ),
    moveMeshOuterCorrectors_
    (
        pimple().dict().lookupOrDefault("moveMeshOuterCorrectors", false)
    ),
    cumulativeContErr_(0),
    adjustTimeStep_
    (
        runTime.controlDict().lookupOrDefault<Switch>("adjustTimeStep", false)
    ),
        maxCo_
    (
        runTime.controlDict().lookupOrDefault<scalar>("maxCo", 1.0)
    ),
    maxDeltaT_
    (
        runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT)
    )
 
{



    UisRequired();
    pisRequired();

        // Reset phi dimensions: compressible
    Info<< "Resetting the dimensions of phi 1 " << endl;
    // Info << "Phi "<< phi()<<endl;
    phi().dimensions().reset(dimVelocity*dimArea*dimDensity);

    // Info << "linearInterpolate(rho_*U()) & mesh().Sf() "<< linearInterpolate(rho_*U()) & mesh().Sf()<<endl;
    Info<< "Resetting the dimensions of phi 2 " << endl;
        // Info << "Phi "<< phi()<<endl;
        
    phi() = linearInterpolate(rho_*U()) & mesh().Sf();
        Info<< "Resetting the dimensions of phi 3 " << endl;

    mesh().setFluxRequired(p().name());
    p() = thermo_.p();

    turbulence_->validate();

    if (mesh().dynamic())
    {
        Info<< "Constructing face velocity Uf\n" << endl;

        rhoUf_.reset
        (
            new surfaceVectorField
            (
                IOobject
                (
                    "rhoUf",
                    runTime.timeName(),
                    mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                fvc::interpolate(rho_*U())
            )
        );


        rhoUf_().oldTime();

    
        if (U().nOldTimes())
        {
            volVectorField* Uold = &U().oldTime();
            volScalarField* Kold = &K_.oldTime();
            *Kold == 0.5*magSqr(*Uold);

            while (Uold->nOldTimes())
            {
                Uold = &Uold->oldTime();
                Kold = &Kold->oldTime();
                *Kold == 0.5*magSqr(*Uold);
            }
        }

    }


    Info<< "exitting pimpleSonicFluid constructor " << endl;

    const fvMesh& mesh = this->mesh();
    const surfaceScalarField& phi = this->phi();
    #include "CourantNo.H"

    // Create temperature field if necessary


 
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void pimpleSonicFluid::setDeltaT(Time& runTime)
{
     Info<< "setDeltaT() frompimpleSonicFluid" <<endl;
    if (adjustTimeStep_)
    {
        Info<< "setDeltaT() from pimpleSonicFluid if (adjustTimeStep_)" <<endl;
        // Calculate the maximum Courant number
        // Careful to use the relative flux in the calculation
        // We have to be careful when we call makeRelative and makeAbsolute
        scalar CoNum = 0.0;
        scalar meanCoNum = 0.0;
        scalar velMag = 0.0;
        fvc::makeRelative(phi(), rho_, U());
        CourantNo(CoNum, meanCoNum, velMag);
        fvc::makeAbsolute(phi(), rho_ , U());

        scalar maxDeltaTFact = maxCo_/(CoNum + SMALL);
        scalar deltaTFact =
            min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

        runTime.setDeltaT
        (
            min
            (
                deltaTFact*runTime.deltaT().value(),
                maxDeltaT_
            )
        );

        Info<< "deltaT = " <<  runTime.deltaT().value() << endl;
    }
}

tmp<vectorField> pimpleSonicFluid::patchViscousForce(const label patchID) const
{
    tmp<vectorField> tvF
    (
        new vectorField(mesh().boundary()[patchID].size(), vector::zero)
    );

    tvF.ref() = 
       (
            mesh().boundary()[patchID].nf()
          & (-turbulence_->devRhoReff()().boundaryField()[patchID])
        );

    return tvF;
}


tmp<scalarField> pimpleSonicFluid::patchPressureForce(const label patchID) const
{
    tmp<scalarField> tpF
    (
        new scalarField(mesh().boundary()[patchID].size(), 0)
    );

    tpF.ref() = p().boundaryField()[patchID];

    return tpF;
}


bool pimpleSonicFluid::evolve()
{
    Info<< "Evolving fluid model: " << this->type() << endl;

    // Take references
    const Time& runTime = fluidModel::runTime();
    dynamicFvMesh& mesh = this->mesh();
    pimpleControl& pimple = this->pimple();
    psiThermo& thermo = thermo_;
    volVectorField& U = this->U();
    volScalarField& p = this->p();
    surfaceScalarField& phi = this->phi();
    autoPtr<surfaceVectorField>& rhoUf = this->rhoUf_;
    volScalarField& K = K_;

    volScalarField& e = e_;

    const volScalarField& psi = psi_;

    //- Fluid density
    volScalarField& rho =  rho_;

    scalar& cumulativeContErr = cumulativeContErr_;
    const bool correctPhi = correctPhi_;
    const bool checkMeshCourantNo = checkMeshCourantNo_;
    const bool moveMeshOuterCorrectors = moveMeshOuterCorrectors_;

    // --- Pressure-velocity PIMPLE corrector loop

    #include "rhoEqn.H"
    Info<< "rho min/max : " << min(rho).value() << " " << max(rho).value()
            << endl;

    while (pimple.loop())
    {
        if (pimple.firstIter() || moveMeshOuterCorrectors)
        {

            // Ideally we would not need a specific FSI mesh update function
            // Hopefully we can remove the need for it soon
            if (fluidModel::fsiMeshUpdate())
            {
                // The FSI interface is in charge of calling mesh.update()
                fluidModel::fsiMeshUpdateChanged();
            }
            else
            {
                // Do any mesh changes
                mesh.controlledUpdate();
            }

            if (mesh.changing())
            {

                if (correctPhi)
                {
                    // Calculate absolute flux
                    // from the mapped surface velocity
                     phi = mesh.Sf() & rhoUf();
                    #include "correctPhi.esi.H"

                    // Make the flux relative to the mesh motion
                    fvc::makeRelative(phi,rho, U);
                }

                if (checkMeshCourantNo)
                {
                    #include "meshCourantNo.H"
                }
            }
        }
    
        #include "UEqn.H"
        #include "EEqn.H"

       // --- Pressure corrector loop
        while (pimple.correct())
        {
            #include "pEqn.H"

        }

        if (pimple.turbCorr())
        {
            turbulence_->correct();
        }
     
    }
    
    rho = thermo.rho();

    return 0;
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels
} // End namespace Foam

// ************************************************************************* //
