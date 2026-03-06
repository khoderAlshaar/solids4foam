{
              
  // TimeState subCycleTimeState = runTime.subCycle(numberSubCycles);

  physDeltaT[1] = runTime.deltaT().value();
  physDeltaT[2] = runTime.deltaT0().value();
  physDeltaT[0] = 1.0;


    // dbnsFlux.meshPhi() = fvc::meshPhi(U);
  

  int iCorr = 0;
  // scalar pseudoTimeStep=0;
  // Low storage Runge-Kutta time integration
  // Pseudo time loop (Correction loop)
   scalar initialRezRho=0;
   scalar rhoResNorm = 0;
   scalar resRho =0;
  // do
  // {
  for (int intIter=0; intIter < numberSubCycles; intIter++)
  {
    // Info << "Pesudo iteration " << intIter+1 << nl;



lduMatrix::debug=0;
blockLduMatrix::debug=0;
    // Update subCycle Time
    // runTime++;

    if (!adjustTimeStep)
    {
        localTimeStep.update(maxCo,adjustTimeStep);
    }

    // Update subCycle Time


    // setpDeltaT(CoDeltaT, mesh, U, thermo, maxCo);
    pseudoTimeStep = dimensionedScalar("pseudoTimeStep", CoDeltaT.dimensions(), min(localTimeStep.CoDeltaT().internalField()));
    
    if (pseudoTimeStep > runTime.deltaT())
    {
      FatalErrorIn("pseudoTime loop)")
            << "physical Time step is smaller than pesudo-time step" << nl
            << " Increase physical-time step or reduce CFL number" << nl
            << exit(FatalError);
    }
    rho.storePrevIter();
    p.storePrevIter();
    rhoU.storePrevIter();
    rhoE.storePrevIter();

    // forAll (beta, i)
    for (int i=0 ; i<2 ; i++)
    {   

    // physDeltaT[0] = beta[i];
      // Solve the approximate Riemann problem for this time step
      dbnsFlux.computeFlux();


      volScalarField rhoFlux = - fvc::div(dbnsFlux.rhoFlux());

      volVectorField rhoUFlux = - fvc::div(dbnsFlux.rhoUFlux());

      volScalarField rhoEFlux = - fvc::div(dbnsFlux.rhoEFlux());

                // // Pseudo time integration
                // solve
                // (
                //     fvm::ddt(rho) == rhoFlux
                // );

                // // time integration
                // solve
                // (
                //     fvm::ddt(rhoU) == rhoUFlux
                // );

                // // time integration
                // solve
                // (
                //     fvm::ddt(rhoE) == rhoEFlux
                // );

      // Pseudo time integration
      fvScalarMatrix rhoEqn
      (
          fvm::ddt(rho) == rhoFlux
      );
      
      fvVectorMatrix rhoUEqn
      (
          fvm::ddt(rhoU) == rhoUFlux
      );      
      fvScalarMatrix rhoEEqn
      (
          fvm::ddt(rhoE) == rhoEFlux
      );

      if ( i == 1)
        {
          resRho =gSum(sqr(rhoEqn.residual())) ; 

          // Info << "rho residual (L2 norm): " << resRho << nl;
          if (intIter == 0) 
          {
            initialRezRho  = resRho;
          }
        }
      rhoEqn.solve();
      rhoUEqn.solve();
      rhoEEqn.solve();

                
      #include "newUpdateFields.H"
      // phi = dbnsFlux.rhoFlux();

        //  rhoResNorm = sqrt(gSum(Foam::sqr(rhoFlux.internalField())) / rhoFlux.size());

      // End outer Loop for Runge - Kutta
    }
    rho = 0.5*(rho.prevIter() + rho);
    rhoU = 0.5*(rhoU.prevIter() + rhoU);
    rhoE = 0.5*(rhoE.prevIter() + rhoE);


    #include "newUpdateFields.H"   


    if (mesh.moving())
    {
        phi = (linearInterpolate(rho*U) & mesh.Sf()) - (fvc::meshPhi(U) *fvc::interpolate(rho)); //?
    }
    else
    {
      phi = linearInterpolate(rho*U) & mesh.Sf();
    }  
      
    bool lastIteration = ( intIter + 1 == numberSubCycles );

          lastIteration = lastIteration || 
              ( resRho < relTol * initialRezRho );

          lastIteration = lastIteration || 
              (resRho < tolerance);

      if (lastIteration) {
          Info << "Pseudo solution  converged in " << intIter+1 << " iterations, " << "residual= "<< resRho<< nl;
          break;
      }
      else if (intIter > numberSubCycles ) {
          Info << "Pseudo solution did not converged in " << intIter+1 << " iterations, " << "residual= "<< resRho<< nl;
          break;
      }


      if (intIter % 100 == 0)
      {
        Info << "iteration No= " << intIter << " residual= "<< resRho << nl;

      }
  }

  // }
  // while
  // (
  //   !converged
  //       (
  //         iCorr,
  //         pseudoTimeStep,
  //         p,
  //         numberSubCycles,
  //         tolerance
  //       )
  // && ++iCorr < numberSubCycles && rhoResNorm < tolerance
  // );
  // Info << " ------------------- solution converged ---- icor = "<<iCorr<<"---------- "<<endl;
    // End dual-time stepping
  // runTime.endSubCycle();
      
    // #include "updateDualTimeSteppingFields.H"

}