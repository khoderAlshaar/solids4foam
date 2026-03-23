{
    // Estimate the bdensity at the boundary using "zeroGradient"
    // rho.boundaryField() = thermo.rho()().boundaryField();
    // Compute U
    U = rhoU/rho;
    U.correctBoundaryConditions();

    volScalarField e = rhoE/rho - 0.5*magSqr(U); 
    e.correctBoundaryConditions();
    // Compute p from rho
    // p = (gamma -1) *rho*(e - q) - gamma * pInf;
    const volScalarField gamma = thermo.Cp()/thermo.Cv();

    p = (gamma - 1.0)*rho*(e - q) - gamma*pInf;
    p.correctBoundaryConditions();
    // // Bound pressure   
    // // boundMinMax(p, pMin, pMax);
    // // h = Cp/Cv*(rhoE/rho - 0.5*magSqr(U));
    // // Info<< "before update h = "<<h<<endl;
    // // h = (rhoE/rho - 0.5*magSqr(U)) + (p/rho);
    // h = (gamma*(p + pInf))/((gamma - 1)*rho);
    h = e + (p/rho) ;
    h.correctBoundaryConditions();


    // T = (p + pInf) /(((gamma-1)/gamma)*cp *rho); 
    // // Correct thermodynamics
    // // Info <<"----------------thermo.correct() start -----------------------------------"<<endl;
    thermo.correct();
    // // Info <<"----------------thermo.correct() end -----------------------------------"<<endl;


    // Update boundary field of rho
    rho.boundaryField() = thermo.rho()().boundaryField();
    // rho.boundaryField() = thermo.psi().boundaryField() * (p.boundaryField() + pInf.boundaryField());

    // Bound density
    // boundMinMax(rho, rhoMin, rhoMax);

    // Update boundary field of rhoU 
    rhoU.boundaryField() = rho.boundaryField()*U.boundaryField();
    
    rhoE.boundaryField() =   rho.boundaryField()*(((p.boundaryField() + (gamma.boundaryField()*pInf.value()))/(p.boundaryField() + pInf.value()) )*thermo.Cv()().boundaryField() * T.boundaryField() + q.value()) +0.5*rho.boundaryField()*magSqr(U.boundaryField());

}
