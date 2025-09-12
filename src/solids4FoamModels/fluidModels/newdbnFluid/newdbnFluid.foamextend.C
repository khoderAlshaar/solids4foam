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

#include "newdbnFluid.H"
// #include "volFields.H"
// #include "fvm.H"
// #include "fvc.H"
// #include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
// #include "findRefCell.H"
// #include "adjustPhi.H"
// #include "wedgeFvPatchFields.H"
// #include "slipFvPatchFields.H"
// #include "extrapolatedFvPatchFields.H"
// #include "EulerDdtScheme.H"
// #include "backwardDdtScheme.H"
// #include "elasticSlipWallVelocityFvPatchVectorField.H"
// #include "elasticWallVelocityFvPatchVectorField.H"
// #include "elasticWallPressureFvPatchScalarField.H"
// #include "movingWallPressureFvPatchScalarField.H"
// #include "flowRateOutletPressureFvPatchScalarField.H"
// #include "thermalRobinFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(newdbnFluid, 0);
addToRunTimeSelectionTable(fluidModel, newdbnFluid, dictionary);

// * * * * * * * * * * * * * * * Private Members * * * * * * * * * * * * * * //

bool newdbnFluid::converged
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
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

newdbnFluid::newdbnFluid
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
    p_ (thermo_.p()),
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
    U_
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh()
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
            // U_,
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
    rho_.oldTime();
    rhoU_.oldTime();
    rhoE_.oldTime();
    rho_.oldTime().oldTime();
    rhoU_.oldTime().oldTime();
    rhoE_.oldTime().oldTime();

    phi().dimensions().reset(dimVelocity*dimArea*dimDensity);

    phi() = (linearInterpolate(rho_*U()) & mesh().Sf()) - (fvc::meshPhi(U()) *linearInterpolate(rho_)); //?

    p() = thermo_.p();

    // turbulence_.validate();

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


tmp<vectorField> newdbnFluid::patchViscousForce(const label patchID) const
{
    tmp<vectorField> tvF
    (
        new vectorField(mesh().boundary()[patchID].size(), vector::zero)
    );

// #ifdef OPENFOAM_NOT_EXTEND
//     tvF.ref() =
// #else
//     tvF() =
// #endif
//         rho_.value()
//        *(
//             mesh().boundary()[patchID].nf()
//           & (-turbulence_->devReff()().boundaryField()[patchID])
//         );
        // rho_.value()*nu_.value()*U_.boundaryField()[patchID].snGrad();

    return tvF;
}


tmp<scalarField> newdbnFluid::patchPressureForce(const label patchID) const
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


bool newdbnFluid::evolve()
{
    Info << "Evolving fluid model: " << this->type() << endl;
    // surfaceScalarField& phi = this->phi();

    // Info << U() <<endl;
    #include "readFieldBounds.H"
    #include "readMultiStage.H"
    
    if (adjustTimeStep_)
    {
        localTimeStep_.update(maxCo_,adjustTimeStep_);
        runTime_.setDeltaT
        (
            min
            (
                min(localTimeStep_.CoDeltaT()).value(),
                maxDeltaT_
            )
        );
        numberSubCycles_ = 1;
    }

    // Info<< "\n physical Time = " << runTime_.value() << endl;
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
        mesh().update();
    }

    #include "mySolveRK2.H"

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels

} // End namespace Foam

// ************************************************************************* //
