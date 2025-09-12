{

                // Compute internal field of U
                U().dimensionedInternalField() =
                    rhoU_.dimensionedInternalField()
                   /rho_.dimensionedInternalField();

                // Update boundary conditions of U
                U().correctBoundaryConditions();

                // store tmp fields in order to prevent computing twice
                const volScalarField kappa(thermo_.Cp()/thermo_.Cv());

                h_ = kappa*(rhoE_/rho_ - 0.5*magSqr(U()));

                // correct boundary conditions of static enthalpy
                h_.correctBoundaryConditions();

                thermo_.correct();
                // compute complete field of p
//                 p = (kappa - 1.0)*rho*e;
                p_ = (1.0 - 1.0/kappa)*rho_*h_;

                // correct boundary conditions of p
                p_.correctBoundaryConditions();

                // Update boundary field of rho
                rho_.boundaryField() = thermo_.rho()().boundaryField();

                // Update boundary field of rhoU
                rhoU_.boundaryField() = rho_.boundaryField()*U().boundaryField();

                // Update boundary field of rhoE
                rhoE_.boundaryField() =
                    rho_.boundaryField()*
                    (
                        0.5*magSqr(U().boundaryField())
                      // + turbulence.k()().boundaryField()
                    )
                  + p_.boundaryField()/(kappa.boundaryField()-1.0);

          p() = p_;
}
