{
    // Compute U
                 // Compute internal field of U
    // U.dimensionedInternalField() =
    //                 rhoU.dimensionedInternalField()
    //                /rho.dimensionedInternalField();
    U = rhoU/rho;

    // Update boundary conditions of U
    U.correctBoundaryConditions();


    const volScalarField kappa(thermo.Cp()/thermo.Cv());


    p = (kappa-1)*(rhoE- 0.5*rho*magSqr(U));
    p.correctBoundaryConditions();

    // h = rhoE/rho - 0.5*magSqr(U) + p/rho;
    // Alternative formulation, Felipe Alves Portela TU Delft
    h = kappa*(rhoE/rho - 0.5*magSqr(U));
    h.correctBoundaryConditions();



    // Bound enthalpy
    // dimensionedScalar CpMin = min(Cp);
    // dimensionedScalar CpMax = max(Cp);

    // dimensionedScalar hMin = CpMin*TMin;
    // dimensionedScalar hMax = CpMax*TMax;

    // boundMinMax(h, hMin, hMax);

    // Correct thermodynamics
    thermo.correct();

    // Bound density
    // boundMinMax(rho, rhoMin, rhoMax);

    // Compute p from rho
    // p = rho*(thermo.Cp() - thermo.Cv())*T;
    // p = (1.0 - 1.0/kappa)*rho*h;
    // p.correctBoundaryConditions();

                // Update boundary field of rho
    rho.boundaryField() = p.boundaryField()*thermo.psi().boundaryField();;
    // rho.boundaryField() = thermo.rho()().boundaryField();
    // rho.boundaryField() = thermo.rho()().boundaryField();

                // Update boundary field of rhoU
    rhoU.boundaryField() = rho.boundaryField()*U.boundaryField();

                // Update boundary field of rhoE
    rhoE.boundaryField() = rho.boundaryField()*( 0.5*magSqr(U.boundaryField()) ) + p.boundaryField()/(kappa.boundaryField()-1.0);
}
