

        // if (adjustTimeStep_)
        // {
        //     localTimeStep_.update(maxCo_, localDt_);
        //     runTime.setDeltaT
        //     (

        //             min(localTimeStep_.CoDeltaT()).value()
        //     );
        //     // conv.pseudoMaxIters() = 1;
        // }

        Info << "deltaT = " << runTime.deltaT().value() << endl;
        Info << "physical Time = " << runTime.value() << endl;

        lduMatrix::debug = 0;


                          // Update primitive boundaries
            p.correctBoundaryConditions();
            U.correctBoundaryConditions();
            h.correctBoundaryConditions();
            // Solve the approximate Riemann problem for this time step
            dbnsFlux.computeFlux();
            
            dimensionedScalar dt = runTime.deltaT();
            
            volScalarField dRho( -dt * ( 
                fvc::ddt(rho) + fvc::div(dbnsFlux.rhoFlux()) 
            ));
            
            volVectorField dRhoU( -dt * (
                fvc::ddt(rhoU)
                + fvc::div(dbnsFlux.rhoUFlux())
                + fvc::div(turbulence_->devRhoReff()) 
            )); 

            volScalarField dRhoE(-dt*( 
                fvc::ddt(rhoE)
                + fvc::div(dbnsFlux.rhoEFlux()) 
                + fvc::div(turbulence_->devRhoReff() & U)
                - fvc::laplacian(turbulence_->alphaEff(), h) 
            ));
           
            rho  += dRho;
            rhoU += dRhoU;
            rhoE += dRhoE;
            
            #include "newUpdateFields.H"

            
            // thermo->correct();
    
            if (mesh.moving())
            {
                phi = dbnsFlux.rhoFlux() - (fvc::meshPhi(U) * fvc::interpolate(rho));
            }
            else
            {
                phi = dbnsFlux.rhoFlux();
            }

        
   
    turbulence_->correct();


