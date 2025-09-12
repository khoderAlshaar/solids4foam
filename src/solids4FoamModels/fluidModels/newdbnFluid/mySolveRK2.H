{
              
  physDeltaT_[1] = runTime_.deltaT().value();
  physDeltaT_[2] = runTime_.deltaT0().value();
  physDeltaT_[0] = 1.0;
  

  int iCorr = 0;
  // Low storage Runge-Kutta time integration
  // Pseudo time loop (Correction loop)
  do
  {


        lduMatrix::debug = 0;
    // Update subCycle Time
    // runTime++;
    // Update subCycle Time
    if (!adjustTimeStep_)
    {
        localTimeStep_.update(maxCo_,adjustTimeStep_);
    }

   
    // setpDeltaT(CoDeltaT, mesh, U, thermo, maxCo);
    pseudoTimeStep_ = dimensionedScalar("pseudoTimeStep", CoDeltaT_.dimensions(), min(localTimeStep_.CoDeltaT().internalField()));
    
    // Info << "pseudoTimeStep_" <<pseudoTimeStep_<<endl;
    // Info << "runTime_.deltaT()" <<runTime_.deltaT()<<endl;
    if (pseudoTimeStep_ > runTime_.deltaT())
    {
      FatalErrorIn("pseudoTime loop)")
            << "physical Time step= " << runTime_.deltaT().value()<< "is smaller than minimum pesudo-time step= "<<pseudoTimeStep_ << nl
            << " Increase physical-time step or reduce CFL number" << nl
            << exit(FatalError);
    }
    
    rho_.storePrevIter();
    rhoU_.storePrevIter();
    rhoE_.storePrevIter();

    // forAll (beta, i)
    for (int i=0 ; i<2 ; i++)
    {   
      
    // physDeltaT[0] = beta[i];
      // Solve the approximate Riemann problem for this time step
      dbnsFlux_.computeFlux();

      // Pseudo time integration
      solve
      (
          fvm::ddt(rho_) == - fvc::div(dbnsFlux_.rhoFlux())
      );
      
      solve
      (
          fvm::ddt(rhoU_) == - fvc::div(dbnsFlux_.rhoUFlux())
      );      
      solve
      (
          fvm::ddt(rhoE_) == - fvc::div(dbnsFlux_.rhoEFlux())
      );


                
      #include "newUpdateFields.H"
      // phi = (linearInterpolate(rho*U) & mesh.Sf()) - (fvc::meshPhi(U) *fvc::interpolate(rho)); //?
      // turbulence->correct();

      // phi = fvc::makeRelative ( phi, rho, U);
        //  rhoResNorm = sqrt(gSum(Foam::sqr(rhoFlux.internalField())) / rhoFlux.size());

      // End outer Loop for Runge - Kutta
    }

    rho_ = 0.5*(rho_.prevIter() + rho_);
    rhoU_ = 0.5*(rhoU_.prevIter() + rhoU_);
    rhoE_ = 0.5*(rhoE_.prevIter() + rhoE_);


    #include "newUpdateFields.H" 
    if (mesh().moving())
    {
        phi() = (linearInterpolate(rho_*U()) & mesh().Sf()) - (fvc::meshPhi(U()) *fvc::interpolate(rho_)); //?
    }
    else
    {
      phi() = linearInterpolate(rho_*U()) & mesh().Sf();
    }

          // turbulence->correct();    

  }
  while
  (
    !converged
        (
          iCorr,
          pseudoTimeStep_,
          rho_,
          numberSubCycles_,
          tolerance_
        )
  && ++iCorr < numberSubCycles_ //&& rhoResNorm < tolerance
  );
  Info << " ------------------- solution converged ---- icor = "<<iCorr<<"---------- "<<endl;
// #       include "updateDualTimeSteppingFields.H"

}