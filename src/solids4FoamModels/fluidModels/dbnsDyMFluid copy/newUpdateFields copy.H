{

                // Compute internal field of U
                U.dimensionedInternalField() =
                    rhoU.dimensionedInternalField()
                   /rho.dimensionedInternalField();

                // Update boundary conditions of U
                U.correctBoundaryConditions();

                // store tmp fields in order to prevent computing twice
                const volScalarField kappa(thermo.Cp()/thermo.Cv());

                h = kappa*(rhoE/rho - 0.5*magSqr(U));

                // correct boundary conditions of static enthalpy
                h.correctBoundaryConditions();

                thermo.correct();
                // compute complete field of p
//                 p = (kappa - 1.0)*rho*e;
                p = (1.0 - 1.0/kappa)*rho*h;

                // correct boundary conditions of p
                p.correctBoundaryConditions();

                // Update boundary field of rho
                rho.boundaryField() = thermo.rho()().boundaryField();

                // Update boundary field of rhoU
                rhoU.boundaryField() = rho.boundaryField()*U.boundaryField();

                // Update boundary field of rhoE
                rhoE.boundaryField() =
                    rho.boundaryField()*
                    (
                        0.5*magSqr(U.boundaryField())
                      // + turbulence.k()().boundaryField()
                    )
                  + p.boundaryField()/(kappa.boundaryField()-1.0);


}
