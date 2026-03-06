// mySolveRK2.H - Pseudo-time RK2 inner loop
// Expects fields and control variables to be visible in scope

physDeltaT[1] = runTime.deltaT().value();
physDeltaT[2] = runTime.deltaT0().value();
physDeltaT[0] = 1.0;

int iCorr = 0;
scalar initialResidual = GREAT; // per-physical-step initial residual

// Pseudo time loop (Correction loop)
do
{
  
    lduMatrix::debug = 0;

    if (!adjustTimeStep)
    {
        localTimeStep.update(maxCo, adjustTimeStep);
    }

    pseudoTimeStep = dimensionedScalar
    (
        "pseudoTimeStep",
        CoDeltaT.dimensions(),
        min(localTimeStep.CoDeltaT().internalField())
    );

    if (pseudoTimeStep > runTime.deltaT())
    {
        FatalErrorIn("pseudoTimeLoop")
            << "physical Time step= " << runTime.deltaT().value()
            << " is smaller than pseudo-time step= " << pseudoTimeStep.value() << nl
            << "Increase physical-time step or reduce CFL number" << nl
            << exit(FatalError);
    }

    rho.storePrevIter();
    rhoU.storePrevIter();
    rhoE.storePrevIter();

    // low-storage RK2 (2 stages)
    for (int i = 0; i < 2; ++i)
    {
        dbnsFlux.computeFlux();

        solve( fvm::ddt(rho) == - fvc::div(dbnsFlux.rhoFlux()) );

        solve
        (
            fvm::ddt(rhoU) == - fvc::div(dbnsFlux.rhoUFlux())
                              - fvc::div(turbulence_->devRhoReff())
        );

        solve
        (
            fvm::ddt(rhoE) == - fvc::div(dbnsFlux.rhoEFlux())
                              - fvc::div(turbulence_->devRhoReff() & U)
                              + fvc::laplacian(turbulence_->alphaEff(), h)
        );

        if( i ==0 )
        {
            #include "newUpdateFields.H"
                    
            if (mesh.moving())
            {
                phi = (linearInterpolate(rho*U) & mesh.Sf()) - (fvc::meshPhi(U) * fvc::interpolate(rho));
            }
            else
            {
                phi = linearInterpolate(rho*U) & mesh.Sf();
            }
        }

    }

    rho  = 0.5*(rho.prevIter()  + rho);
    rhoU = 0.5*(rhoU.prevIter() + rhoU);
    rhoE = 0.5*(rhoE.prevIter() + rhoE);

    #include "newUpdateFields.H"


    if (mesh.moving())
    {
        phi = (linearInterpolate(rho*U) & mesh.Sf()) - (fvc::meshPhi(U) * fvc::interpolate(rho));
    }
    else
    {
        phi = linearInterpolate(rho*U) & mesh.Sf();
    }


}
while (!conv.pseudoConverged(iCorr, pseudoTimeStep, rho, rhoU, rhoE) && ++iCorr < conv.pseudoMaxIters());

turbulence_->correct();

// Info << " ------------------- pseudo-time converged, iCorr = " << iCorr << " ---------- " << nl;
