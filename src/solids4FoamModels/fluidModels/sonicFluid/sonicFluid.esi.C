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

#include "sonicFluid.H"
#include "addToRunTimeSelectionTable.H"
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

defineTypeNameAndDebug(sonicFluid, 0);
addToRunTimeSelectionTable(fluidModel, sonicFluid, dictionary);


// * * * * * * * * * * * * * * * Private Members * * * * * * * * * * * * * * //

void sonicFluid::compressibleContinuityErrs()
{
    dimensionedScalar totalMass = fvc::domainIntegrate(rho_);

    scalar sumLocalContErr =
        (fvc::domainIntegrate(mag(rho_ - thermo_.rho()))/totalMass).value();

    scalar globalContErr =
        (fvc::domainIntegrate(rho_ - thermo_.rho())/totalMass).value();

    cumulativeContErr_ += globalContErr;

    Info<< "time step continuity errors : sum local = " << sumLocalContErr
        << ", global = " << globalContErr
        << ", cumulative = " << cumulativeContErr_
        << endl;
}


void sonicFluid::solveRhoEqn()
{
    fvScalarMatrix rhoEqn
    (
        fvm::ddt(rho_)
      + fvc::div(phi())
     ==
        options()(rho_)
    );

    options().constrain(rhoEqn);

    rhoEqn.solve();

    options().correct(rho_);
}


void sonicFluid::CorrectFlux()
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
                fvc::div(phi())
        )
    );

    // Store momentum to set rhoUf for introduced faces.
    autoPtr<volVectorField> rhoU;

    // Calculate absolute flux
    // from the mapped surface velocity
    phi() = mesh().Sf() & rhoUf_();

    // Define volScalarField to hold phi
    // to pass to compressible CorrectPhi function
    // const volScalarField psi
    // (
    //      IOobject
    //      (
    //             "psi",
    //             runTime().timeName(),
    //             mesh(),
    //             IOobject::NO_READ,
    //             IOobject::NO_WRITE
    //      ),
    //      mesh(),
    //      psi_
    // );

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

void sonicFluid::compressibleCourantNo()
{
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;
    scalar velMag = 0.0;

    // HR 26.06.18: A parallel run has at least two cells and therefore at least
    // one internal face in the global mesh. It may be a processor boundary, but
    // this is captured by max(mag(phi)).
    // Old formulation hangs on parallel cases where one partition is degenerated
    // to a single cell.
    if (mesh().nInternalFaces() || Pstream::parRun())
    {
        surfaceScalarField phiOverRho = mag(phi())/fvc::interpolate(rho_);

        surfaceScalarField SfUfbyDelta =
            mesh().surfaceInterpolation::deltaCoeffs()*phiOverRho;

        CoNum = max(SfUfbyDelta/mesh().magSf()).value()*runTime().deltaT().value();

        meanCoNum = (sum(SfUfbyDelta)/sum(mesh().magSf())).value()*
            runTime().deltaT().value();

        velMag = max(phiOverRho/mesh().magSf()).value();
    }


    Info<< "Courant Number mean: " << meanCoNum
        << " max: " << CoNum
        << " velocity magnitude: " << velMag
        << endl;
}


void sonicFluid::solvePEqn(const fvVectorMatrix& UEqn)
{

    rho_ = thermo_.rho();
    volScalarField rAU(1.0/UEqn.A());   

    const surfaceScalarField rhorAUf
    (
        "rhorAUf",
        fvc::interpolate(rho_*rAU)
    );

    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U(), p())); //!new

    surfaceScalarField phid
    (
        "phid",
        fvc::interpolate(psi_)
        *(
            fvc::flux(HbyA)
            // + rhorAUf*fvc::ddtCorr(rho_, U(), rhoUf_())/fvc::interpolate(rho_)
            + rhorAUf*fvc::ddtCorr(rho_, U(), phi())/fvc::interpolate(rho_)
            )
    );
    

    // Make flux relative to mesh motion
    // fvc::makeRelative(phid, psi_, U());

// Non-orthogonal pressure corrector loop
    while (pimple().correctNonOrthogonal())
    {
        // Pressure equation
        //! essential to account for ddt(rho0) term for mesh motion)
        fvScalarMatrix pEqn
        (
              fvm::ddt(psi_, p())
            + fvm::div(phid, p())
            - fvm::laplacian(rhorAUf, p())
           ==
                options()(psi_, p(), rho_.name())
        );

        pEqn.solve();

        if (pimple().finalNonOrthogonalIter())
        {
            phi() += pEqn.flux();
        }
    }

    solveRhoEqn();

    // State equation errors
    compressibleContinuityErrs();

    // Correct velocity

    U() = HbyA - rAU*fvc::grad(p());
    U().correctBoundaryConditions();
    options().correct(U());
    K_ = 0.5*magSqr(U());

    // if (mesh().dynamic())
    // {
    //     rhoUf_() = fvc::interpolate(rho_*U());
    //     surfaceVectorField n(mesh().Sf()/mesh().magSf());
    //     rhoUf_() += n*(fvc::absolute(phi(), rho_, U())/mesh().magSf() - (n & rhoUf_()));

    // }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

sonicFluid::sonicFluid
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
    
    MRF(mesh()),


    cumulativeContErr_(0)
 { 

    // thermo_.validate(args.executable(), "e");

    p() = thermo_.p();

    UisRequired();
    // pisRequired();

    // Reset phi dimensions: compressible
    Info<< "Resetting the dimensions of phi" << endl;

    phi().dimensions().reset(dimVelocity*dimArea*dimDensity);

    phi() = linearInterpolate(rho_*U()) & mesh().Sf();

    mesh().setFluxRequired(p().name());

    if (mesh().dynamic())
    {
        Info<< "Constructing face momentum rhoUf" << endl;

        rhoUf_.set
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
    }

        turbulence_ = compressible::turbulenceModel::New
    (
        rho_,
        U(),
        phi(),
        thermo_
    );

    
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

    turbulence_->validate();

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<vectorField> sonicFluid::patchViscousForce
(
    const label patchID
) const
{
    tmp<vectorField> tvF
    (
        new vectorField(mesh().boundary()[patchID].size(), vector::zero)
    );

    tvF.ref() = mu_*U().boundaryField()[patchID].snGrad();

    return tvF;
}


tmp<scalarField> sonicFluid::patchPressureForce
(
    const label patchID
) const
{
    tmp<scalarField> tpF
    (
        new scalarField(mesh().boundary()[patchID].size(), 0)
    );

    // Pressure here is already in Pa
    tpF.ref() = p().boundaryField()[patchID];

    return tpF;
}


bool sonicFluid::evolve()
{
    Info<< "Evolving fluid model: " << this->type() << endl;

    // dynamicFvMesh& mesh = this->mesh();

    // bool meshChanged = false;

    // Info <<"runTime().deltaT().value()" <<runTime().deltaT().value()<<endl;

//     if (fluidModel::fsiMeshUpdate())
//     {
//         // The FSI interface is in charge of calling mesh.update()
//         meshChanged = fluidModel::fsiMeshUpdateChanged();
//     }
//     else
//     {
//         meshChanged = mesh.update();
//         reduce(meshChanged, orOp<bool>());
//     }

//     if (meshChanged)
//     {
//         const Time& runTime = fluidModel::runTime();
// #       include "volContinuity.H"
//     }

//     bool correctPhi
//     (
//         pimple().dict().lookupOrDefault("correctPhi", false)
//     );

//     if (correctPhi && meshChanged)
//     {
//         CorrectFlux();
//     }

//     // Make the fluxes relative to the mesh motion
//     fvc::makeRelative(phi(), rho_, U());

    // Calculate CourantNo
    compressibleCourantNo();
        
        // solveRhoEqn();
    
    solve(fvm::ddt(rho_) + fvc::div(phi()));
   
    // Pressure-velocity corrector
    while (pimple().loop())
    {

            #include "UEqn.H"
            #include "EEqn.H"
        // --- Pressure corrector loop
        while (pimple().correct())
        {
            #include "pEqn.H"

        }

        if (pimple().turbCorr())
        {
            turbulence_->correct();
        }
     
        // tUEqn.clear();

    }

    rho_ = thermo_.rho();
    // Make the fluxes absolute to the mesh motion
    // fvc::makeAbsolute(phi(), rho_, U());

    // Print variables for inspection
    // Density variation
    // scalar deltaRho = max(rho_).value() - min(rho_).value();
    // scalar refDeltaRho = 0.01*rho0_.value();

    // Info variable values
    Info<< nl << "Density: min " << min(rho_).value()
        << " max " << max(rho_).value() << endl;

    // Info<< "Density variation: " << deltaRho
    //     << " ref: " << refDeltaRho << endl;

    Info<< "Pressure: min " << min(p()).value()
        << " max " << max(p()).value() << endl;

    Info<< "Velocity: min " << min(mag(U())).value()
        << " max " << max(mag(U())).value() << nl << nl;

    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels
} // End namespace Foam

// ************************************************************************* //
