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

#include "comDbnsDyMFluid.H"
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

namespace compressibleFluidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(comDbnsDyMFluid, 0);
addToRunTimeSelectionTable(compressibleFluidModel, comDbnsDyMFluid, dictionary);

// * * * * * * * * * * * * * * * Private Members * * * * * * * * * * * * * * //




// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

comDbnsDyMFluid::comDbnsDyMFluid
(
    Time& runTime,
    const word& region
)
:
    compressibleFluidModel(typeName, runTime, region),
    runTime_(runTime),
    pThermo_
    (
        basicPsiThermo::New(mesh())
    ),
     thermo_(pThermo_()),
    h_ (thermo_.h()),
    T_ (thermo_.T()),
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
    ),
    
    // dbnsFluxPtr_ 
    // (
    //     // basicNumericFlux::New
    //     // (
            // p(),
            // U(),
            // T_,
            // thermo_
    //     // )
    // ),
    
    // dbnsFlux_ (dbnsFluxPtr_()),
    
    dbnsFluxPtr_ 
    (
        new numericFlux
        (
            p(),
            U(),
            T_,
            thermo_
        )
    ),
    
    dbnsFlux_ (dbnsFluxPtr_()),


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
    
    CoDeltaT_
    (
        IOobject
        (
            "CoDeltaT_",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("CoDeltaT_", dimTime, 0.1)
    ),
    pseudoTimeStep_ ("pseudoTimeStep", dimTime, 0.0),
    numberSubCycles_
    (
        fluidProperties().lookupOrDefault<label>("numberSubCycles", 1)
    ),
    tolerance_
    (
        fluidProperties().lookupOrDefault<scalar>("tolerance", 1e-6)
    ),
    relTol_
    (
        fluidProperties().lookupOrDefault<scalar>("relTol", 1e-6)
    ),

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

    // turbulence_.validate();

    // if (mesh().dynamic())
    // {
    //     Info<< "Constructing face velocity Uf\n" << endl;

    //     rhoUf_.reset
    //     (
    //         new surfaceVectorField
    //         (
    //             IOobject
    //             (
    //                 "rhoUf",
    //                 runTime.timeName(),
    //                 mesh(),
    //                 IOobject::READ_IF_PRESENT,
    //                 IOobject::AUTO_WRITE
    //             ),
    //             fvc::interpolate(rho_*U())
    //         )
    //     );


    //     rhoUf_().oldTime();

    
    //     // if (U().nOldTimes())
    //     // {
    //     //     volVectorField* Uold = &U().oldTime();
    //     //     volScalarField* Kold = &K_.oldTime();
    //     //     *Kold == 0.5*magSqr(*Uold);

    //     //     while (Uold->nOldTimes())
    //     //     {
    //     //         Uold = &Uold->oldTime();
    //     //         Kold = &Kold->oldTime();
    //     //         *Kold == 0.5*magSqr(*Uold);
    //     //     }
    //     // }

    // }


    // Info<< "exitting comDbnsDyMFluid constructor " << endl;

    // const fvMesh& mesh = this->mesh();
    // const surfaceScalarField& phi = this->phi();
    // #include "CourantNo.H"

    // Create temperature field if necessary


 
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool comDbnsDyMFluid::converged
(
    const int iCorr,
    const dimensionedScalar pDeltaT,
    // const volVectorField& vf
    const volScalarField& vf,
    const label numberSubCycles,
    scalar tolerance   
)
{
    // We will check three residuals:
    // - relative linear momentum residual

    bool converged = false;

    // Calculate residual based on the relative change of vf
    scalar denom = 0.0;

    // Denom is linear momentum increment
        denom = gMax
        (
#ifdef OPENFOAM_NOT_EXTEND
            DimensionedField<scalar, volMesh>
#else
            Field<scalar>
#endif
            (
                mag(vf.internalField() - vf.oldTime().internalField())
                // mag(vf.internalField())
            )
        );

    if (denom < SMALL)
    {
        denom =
        max
        (
            gMax
            (
#ifdef OPENFOAM_NOT_EXTEND
                DimensionedField<scalar, volMesh>(mag(vf.internalField()))
#else
                mag(vf.internalField())
#endif
            ),
            SMALL
        );
    }

    const scalar residualvf =
        gMax
        (
#ifdef OPENFOAM_NOT_EXTEND
            DimensionedField<scalar, volMesh>
            (
                mag(vf.internalField() - vf.prevIter().internalField())
            )
#else
            mag(vf.internalField() - vf.prevIter().internalField())
#endif
        )/denom;

    // If one of the residuals has converged to an order of magnitude
    // less than the tolerance then consider the solution converged
    // force at least 1 outer iteration
   if (residualvf < tolerance)
    {
        Info<< "    Converged" << endl;
        converged = true;
    }

    // Print residual information
    if (iCorr == 0)
    {
        Info<< "    Corr, res, pDeltaT" << endl;
    }
    else if (iCorr % 10 == 0 || converged || iCorr >= numberSubCycles)
    {
        Info<< "    " << iCorr
            << ", " << residualvf
            << ", " << pDeltaT.value() << endl;

        if (iCorr >= numberSubCycles)
        {
            Warning
                << "Max iterations reached within the momentum loop"
                << endl;
            converged = true;
        }
    }

    return converged;
 }
void comDbnsDyMFluid::setDeltaT(Time& runTime)
{
     Info<< "comDbnsDyMFluid::setDeltaT(Time& runTime) is off" <<endl;
       
        // if (adjustTimeStep)
        // {
        //     localTimeStep_.update(maxCo_,adjustTimeStep_);
        //     runTime.setDeltaT
        //     (
        //         min
        //         (
        //             min(localTimeStep_.CoDeltaT()).value(),
        //             maxDeltaT
        //         )
        //     );
        //     numberSubCycles_ = 1;
        // }

        // Info<< "\n physical Time = " << runTime.value() << endl;
    
}

tmp<vectorField> comDbnsDyMFluid::patchViscousForce(const label patchID) const
{
    tmp<vectorField> tvF
    (
        new vectorField(mesh().boundary()[patchID].size(), vector::zero)
    );

    // tvF.ref() = 
    //    (
    //         mesh().boundary()[patchID].nf()
    //       & (-turbulence_->devRhoReff()().boundaryField()[patchID])
    //     );

    return tvF;
}

tmp<scalarField> comDbnsDyMFluid::patchPressureForce(const label patchID) const
{
    tmp<scalarField> tpF
    (
        new scalarField(mesh().boundary()[patchID].size(), 0)
    );

#ifdef OPENFOAM_NOT_EXTEND
    tpF.ref() =
#else
    tpF() =
#endif
        p().boundaryField()[patchID];

    return tpF;
}

bool comDbnsDyMFluid::evolve()
{
    Info<< "Evolving fluid model: " << this->type() << endl;

    // // Take references
    Time& runTime = runTime_;
    dynamicFvMesh& mesh = this->mesh();
    basicPsiThermo& thermo = thermo_;
    // basicNumericFlux& dbnsFlux = dbnsFlux_;
    numericFlux& dbnsFlux = dbnsFlux_;
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
    volScalarField& CoDeltaT = CoDeltaT_;



    label& numberSubCycles = numberSubCycles_;
    scalar& tolerance = tolerance_;
    scalar& relTol = relTol_;

    const Switch& adjustTimeStep = adjustTimeStep_;

            //- For adjustable time-step, this is the maximum Courant number
    const scalar& maxCo = maxCo_;

    //- For adjustable time-step, this is the maximum time-stpe
    const scalar& maxDeltaT = maxDeltaT_;
    surfaceScalarField& phi = this->phi();

    #include "readFieldBounds.H"
    #include "readMultiStage.H"

    if (adjustTimeStep)
    {
        localTimeStep.update(maxCo,adjustTimeStep);
        runTime.setDeltaT
        (
            min
            (
                min(localTimeStep.CoDeltaT()).value(),
                maxDeltaT
            )
        );
        numberSubCycles = 1;
    }

    Info<< "\n physical Time = " << runTime.value() << endl;
    // Ideally we would not need a specific FSI mesh update function
    // Hopefully we can remove the need for it soon
    if (compressibleFluidModel::fsiMeshUpdate())
    {
        // The FSI interface is in charge of calling mesh.update()
        compressibleFluidModel::fsiMeshUpdateChanged();
    }
    else
    {
        // Do any mesh changes
        mesh.update();
    }

    // #include "mySolve.H"
    #include "mySolveRK2New.H"
    // #include "mySolveRK2.H"

    return 0;
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressibleFluidModels
} // End namespace Foam

// ************************************************************************* //
