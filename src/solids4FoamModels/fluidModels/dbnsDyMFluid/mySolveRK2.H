{
              
  physDeltaT[1] = runTime.deltaT().value();
  physDeltaT[2] = runTime.deltaT0().value();
  physDeltaT[0] = 1.0;
  

  int iCorr = 0;
  // Low storage Runge-Kutta time integration
  // Pseudo time loop (Correction loop)
  do
  {


        lduMatrix::debug = 0;
    // Update subCycle Time
    // runTime++;
    // Update subCycle Time
    if (!adjustTimeStep)
    {
        localTimeStep.update(maxCo,adjustTimeStep);
    }

   
    // setpDeltaT(CoDeltaT, mesh, U, thermo, maxCo);
    pseudoTimeStep = dimensionedScalar("pseudoTimeStep", CoDeltaT.dimensions(), min(localTimeStep.CoDeltaT().internalField()));
    if (pseudoTimeStep > runTime.deltaT())
    {
      FatalErrorIn("pseudoTime loop)")
            << "physical Time step= "<<runTime.deltaT().value()<<" is smaller than pesudo-time step= "<< pseudoTimeStep.value()<< nl
            << " Increase physical-time step or reduce CFL number" << nl
            << exit(FatalError);
    }
    
    rho.storePrevIter();
    rhoU.storePrevIter();
    rhoE.storePrevIter();

    // forAll (beta, i)
    for (int i=0 ; i<2 ; i++)
    {   
      
    // physDeltaT[0] = beta[i];
      // Solve the approximate Riemann problem for this time step
      dbnsFlux.computeFlux();


      // Pseudo time integration
      solve
      (
          fvm::ddt(rho) == - fvc::div(dbnsFlux.rhoFlux())
      );
      
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


                
      #include "newUpdateFields.H"

      // End outer Loop for Runge - Kutta
    }

    rho = 0.5*(rho.prevIter() + rho);
    rhoU = 0.5*(rhoU.prevIter() + rhoU);
    rhoE = 0.5*(rhoE.prevIter() + rhoE);


      #include "newUpdateFields.H"
      turbulence_->correct();
    if (mesh.moving())
    {
        phi = (linearInterpolate(rho*U) & mesh.Sf()) - (fvc::meshPhi(U) *fvc::interpolate(rho)); //?
    }
    else
    {
      phi = linearInterpolate(rho*U) & mesh.Sf();
    }

          // turbulence->correct();    

  }
while (!conv.pseudoConverged(iCorr, pseudoTimeStep, rho, rhoU, rhoE) && ++iCorr < conv.pseudoMaxIters());
turbulence_->correct();

}