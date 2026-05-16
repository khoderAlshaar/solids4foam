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

#include "dbnsDyMFluid.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "fvm.H"
#include "elasticSlipWallVelocityFvPatchVectorField.H"
#include "elasticWallVelocityFvPatchVectorField.H"
#include "elasticWallPressureFvPatchScalarField.H"
#include "movingWallPressureFvPatchScalarField.H"
#include "EulerDdtScheme.H"
#include "backwardDdtScheme.H"
#include "thermalRobinFvPatchScalarField.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(dbnsDyMFluid, 0);
addToRunTimeSelectionTable(fluidModel, dbnsDyMFluid, dictionary);

// * * * * * * * * * * * * * * * Private Members * * * * * * * * * * * * * * //




// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

dbnsDyMFluid::dbnsDyMFluid
(
    Time& runTime,
    const word& region
)
:
    fluidModel(typeName, runTime, region),
    runTime_(runTime),
    pThermo_
    (
        basicPsiThermo::New(mesh())
    ),
     thermo_(pThermo_()),
    h_ (thermo_.h()),
    T_ (thermo_.T()),
    
    thermoDict
    (
        IOobject
        (
            "thermophysicalProperties",
           runTime.constant(),
            mesh(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    pInf
    (
        "pInf",
        dimPressure,
        thermoDict.subDict("stiffenedGasCoeffs").lookupOrDefault<scalar>("pInf", 0.0)
    ),
    q
    (
        "q",
         dimEnergy/dimMass,
        thermoDict.subDict("stiffenedGasCoeffs").lookupOrDefault<scalar>("q", 0.0)
    ),
    pRef
    (
        "pRef",
         dimPressure,
        thermoDict.subDict("stiffenedGasCoeffs").lookupOrDefault<scalar>("pRef", 0.0)
    ),
    // gamma
    // (
    //     "gamma",
    //     dimless,
    //     thermoDict.subDict("stiffenedGasCoeffs").lookupOrDefault<scalar>("gamma",0.0)
    // ),
    
    rho_
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo_.rho()
    ),  

    rhoU_
    (
        IOobject
        (
            "rhoU",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho_*U()
    ),
    rhoE_
    (
        IOobject
        (
            "rhoE",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho_*(h_ + 0.5*magSqr(U())) - p()
        // rho_*(((p() + (thermo_.Cp()/thermo_.Cv())*pInf)/(p() + pInf) )*thermo_.Cv() * T_ + q) + 0.5*rho_*magSqr(U())
    ),
    
    dbnsFluxPtr_ 
    (
        basicNumericFlux::New
        (
            p(),
            U(),
            T_,
            thermo_
        )
    ),
    
    dbnsFlux_ (dbnsFluxPtr_()),
    
    // dbnsFluxPtr_ 
    // (
    //     new numericFlux
    //     (
    //         p(),
    //         U(),
    //         T_,
    //         thermo_
    //     )
    // ),
    
    // dbnsFlux_ (dbnsFluxPtr_()),


    turbulence_
    (
        compressible::myTurbulenceModel::New
        (
            rho_,
            U(),
            phi(),
            thermo_
        )
    ),
    physDeltaT_
    (
        IOobject
        (
            "physDeltaT",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        3
    ),
    localTimeStep_(U(), thermo_, turbulence_()),
    
    // CoDeltaT_
    // (
    //     IOobject
    //     (
    //         "CoDeltaT_",
    //         runTime.timeName(),
    //         mesh(),
    //         IOobject::NO_READ,
    //         IOobject::NO_WRITE
    //     ),
    //     mesh(),
    //     dimensionedScalar("CoDeltaT_", dimTime, 0.1)
    // ),
    pseudoTimeStep_ ("pseudoTimeStep", dimTime, 0.0),
    // numberSubCycles_
    // (
    //     fluidProperties().lookupOrDefault<label>("numberSubCycles", 1)
    // ),
    // tolerance_
    // (
    //     fluidProperties().lookupOrDefault<scalar>("tolerance", 1e-6)
    // ),
    // relTol_
    // (
    //     fluidProperties().lookupOrDefault<scalar>("relTol", 1e-6)
    // ),
    conv(runTime, mesh()),

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
    ),
    
    localDt_
    (  
        runTime.controlDict().lookupOrDefault("localDt", false)
    ),
    pPrim_
    (
        IOobject
        (
            "pPrim",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p()-pRef
    )
    
{


    conv.read(mesh().solutionDict());
    
    // UisRequired();
    // pisRequired();

    rho_.oldTime();
    rhoU_.oldTime();
    rhoE_.oldTime();
    rho_.oldTime().oldTime();
    rhoU_.oldTime().oldTime();
    rhoE_.oldTime().oldTime();

    //     // Reset phi dimensions: compressible
    // Info<< "Resetting the dimensions of phi 1 " << endl;
    // // Info << "Phi "<< phi()<<endl;
    phi().dimensions().reset(dimVelocity*dimArea*dimDensity);

    // // Info << "linearInterpolate(rho_*U()) & mesh().Sf() "<< linearInterpolate(rho_*U()) & mesh().Sf()<<endl;
    // Info<< "Resetting the dimensions of phi 2 " << endl;
    //     // Info << "Phi "<< phi()<<endl;
        
    // phi() = dbnsFlux_.rhoFlux();     //! this is incorrect due to the use of RK stages
    phi() = (linearInterpolate(rho_*U()) & mesh().Sf()) - (fvc::meshPhi(U()) *linearInterpolate(rho_)); //?
    //     Info<< "Resetting the dimensions of phi 3 " << endl;

    // mesh().setFluxRequired(p().name());
    p() = thermo_.p();


}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void dbnsDyMFluid::setDeltaT(Time& runTime)
{
        // if (adjustTimeStep_)
        // {
        //     localTimeStep_.update(maxCo_, localDt_);
        //     runTime.setDeltaT
        //     (

        //             min(localTimeStep_.CoDeltaT()).value()
        //     );
        //     // conv.pseudoMaxIters() = 1;
        // }
        if (adjustTimeStep_)
        {
        surfaceScalarField amaxSf("amaxSf", 
        mag(fvc::interpolate(U()) & mesh().Sf()) +
        mesh().magSf() * fvc::interpolate(sqrt(thermo_.Cp()/thermo_.Cv()/thermo_.psi())));
        

        #include "compressibleCFLNo.H"
        // #include "setDeltaT.H"
  
            scalar maxDeltaTFact = maxCo_/(CoNum + SMALL);
            scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

            runTime.setDeltaT
            (
                // Foam::max
                // (
                //     minDeltaT_,
                //     Foam::min
                //     (
                        deltaTFact*runTime.deltaT().value()//,
                //         maxDeltaT_
                //     )
                // )
            );

            Info<< "deltaT = " <<  runTime.deltaT().value() << endl;
        }
}



tmp<vectorField> dbnsDyMFluid::patchViscousForce(const label patchID) const
{
    tmp<vectorField> tvF
    (
        new vectorField(mesh().boundary()[patchID].size(), vector::zero)
    );

#ifdef OPENFOAM_NOT_EXTEND
    tvF.ref() =
#else
    tvF() =
#endif
       (
            mesh().boundary()[patchID].nf()
          & (-turbulence_->devRhoReff()().boundaryField()[patchID])
        );
    return tvF;
}


tmp<scalarField> dbnsDyMFluid::patchPressureForce(const label patchID) const
{

    // scalar pRef =1e5;
    tmp<scalarField> tpF
    (
        new scalarField(mesh().boundary()[patchID].size(), 0)
    );

#ifdef OPENFOAM_NOT_EXTEND
    tpF.ref() =
#else
    tpF() =
#endif
        // p().boundaryField()[patchID]- pRef.value();
        p().boundaryField()[patchID];

    return tpF;
}

bool dbnsDyMFluid::evolve()
{
    Info<< "\n Evolving fluid model: " << this->type() << endl;

    // // Take references
    Time& runTime = runTime_;
    dynamicFvMesh& mesh = this->mesh();
    basicPsiThermo& thermo = thermo_;
    basicNumericFlux& dbnsFlux = dbnsFlux_;
    // numericFlux& dbnsFlux = dbnsFlux_;
    volVectorField& U = this->U();
    volScalarField& p = this->p();
    volScalarField& h = this->h_;
    const volScalarField& T = this->T_;
    volScalarField& rho= rho_;
    volVectorField& rhoU = rhoU_;
    volScalarField& rhoE = rhoE_;
    localTimeStep& localTimeStep = localTimeStep_;
    IOField<scalar>& physDeltaT = physDeltaT_;
    dimensionedScalar& pseudoTimeStep = pseudoTimeStep_;
    // volScalarField& CoDeltaT = CoDeltaT_;


    const Switch& adjustTimeStep = adjustTimeStep_;

            //- For adjustable time-step, this is the maximum Courant number
     scalar& maxCo = maxCo_;

    //- For adjustable time-step, this is the maximum time-stpe
    const scalar& maxDeltaT = maxDeltaT_;
    surfaceScalarField& phi = this->phi();

    lduMatrix::debug=0;
    blockLduMatrix::debug=0;

    maxCo = runTime.controlDict().lookupOrDefault<scalar>("maxCo", 1.0);    
 
    #include "readFieldBounds.H"   
    #include "readMultiStage.H"
    conv.read(mesh.solutionDict());

    // if (adjustTimeStep)
    // {
    //     localTimeStep.update(maxCo,adjustTimeStep);
    //     runTime.setDeltaT
    //     (
    //         min
    //         (
    //             min(localTimeStep.CoDeltaT()).value(),
    //             maxDeltaT
    //         )
    //     );
    //     conv.pseudoMaxIters() = 1;
    // }

    // Info<< "\n physical Time = " << runTime.value() << endl;
    // Ideally we would not need a specific FSI mesh update function
    // Hopefully we can remove the need for it soon
            // Info<< "HEllo from 0"<<endl;

    bool meshChanged = false;
    if (fluidModel::fsiMeshUpdate())
    {
        // Info<< "HEllo from 1"<<endl;
        // Info<< "fluidModel::fsiMeshUpdate() "<< fluidModel::fsiMeshUpdate()<<endl;;
        // The FSI interface is in charge of calling mesh.update()
        meshChanged = fluidModel::fsiMeshUpdateChanged(); //! the porblem is that this one stays 0
        // Info<< "meshChanged = fluidModel::fsiMeshUpdateChanged();: "<< meshChanged<<endl;
    }
    else
    {
        //         Info<< "HEllo from 2"<<endl;

        // Info<< "else "<< fluidModel::fsiMeshUpdate();
        meshChanged = mesh.update();
        // Info<< "meshChanged = mesh.update(): "<< meshChanged<<endl;
        reduce(meshChanged, orOp<bool>());
            // Info<< "reduce(meshChanged, orOp<bool>()): "<< meshChanged<<endl;

    }
    // Info<< "meshChanged: "<< meshChanged<<endl;
    if (meshChanged)
    {
                // Info<< "HEllo from 3"<<endl;

        // Info<< "mesh.moving(): "<<mesh.moving()<<endl;
        const Time& runTime = fluidModel::runTime();
#       include "volContinuity.H"
    }
        // Info<< "HEllo from 4"<<endl;

 // Pseudo-time RK2 inner loop (updates rho, rhoU, rhoE, phi, etc.)
        if(conv.solverType() == "MSMS" ) // Multi-Stage-Multi-Step
        {
            #include "mySolveRK2.H"
        }
        else if (conv.solverType() == "MSSS" ) // Multi-Stage-Single-Step
        {
            Info << "solver Type: MSSS -  Multi-Stage-Single-Step"<<endl;
            #include "mySolveMSSS.H"
        }
        else if (conv.solverType() == "SSSS" ) // Single-Stage-Single-Step
        {
            Info << "solver Type: SSSS -  Single-Stage-Single-Step"<<endl;
            #include "mySolveSSSS.H"
        }
        // else if (conv.solverType() == "MSMS-D" ) // Single-Stage-Single-Step
        // {
        //     #include "solveFluid.H"
        // }
        else
        {
            Info << "please set the solver Type"<<endl;
            // break;
        }
        pPrim_ = p - pRef;

    //! After inner loop, check global/physical convergence using oldTime() fields
    conv.physicalConverged(rho, rhoU, rhoE);
    // if (conv.physicalConverged(rho, rhoU, rhoE))
    // {
    //     runTime.write();
    //     // break;
    // }
    
    
    return true;
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels
} // End namespace Foam

// ************************************************************************* //
