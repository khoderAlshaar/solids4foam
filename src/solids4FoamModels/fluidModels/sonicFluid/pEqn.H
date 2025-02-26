rho_ = thermo_.rho();

volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho_*rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U(), p()));
surfaceScalarField phid
(
    "phid",
    fvc::interpolate(psi_)
   *(
        fvc::flux(HbyA)
      + MRF.zeroFilter(rhorAUf*fvc::ddtCorr(rho_, U(), phi())/fvc::interpolate(rho_))
    )
);

MRF.makeRelative(fvc::interpolate(psi_), phid);

// Non-orthogonal pressure corrector loop
while (pimple().correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
        fvm::ddt(psi_, p())
      + fvm::div(phid, p())
      - fvm::laplacian(rhorAUf, p())
     ==
        options()(psi_, p(), rho_.name())
    );

    pEqn.solve();

    if (pimple().finalNonOrthogonalIter())
    {
        phi() = pEqn.flux();
    }
}

solve(fvm::ddt(rho_) + fvc::div(phi()));

compressibleContinuityErrs();

U() = HbyA - rAU*fvc::grad(p());
U().correctBoundaryConditions();
options().correct(U());
K_ = 0.5*magSqr(U());
