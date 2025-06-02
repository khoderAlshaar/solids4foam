rho = thermo_.rho();

volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(),U,p));

    surfaceScalarField phid
    (
        "phid",
        fvc::interpolate(psi)
    *(
            fvc::flux(HbyA)
        + rhorAUf*fvc::ddtCorr(rho,U, phi)/fvc::interpolate(rho)
        )
    );

fvc::makeRelative(phid, psi,U);


// Non-orthogonal pressure corrector loop
while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
        fvm::ddt(psi,p)
      + fvm::div(phid,p)
      - fvm::laplacian(rhorAUf,p)
     ==
        options()(psi,p, rho.name())
    );

    pEqn.solve();

    if (pimple.finalNonOrthogonalIter())
    {
        phi = pEqn.flux();
    }
}

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

U = HbyA - rAU*fvc::grad(p);
U.correctBoundaryConditions();
options().correct(U);
K = 0.5*magSqr(U);

    if (mesh.dynamic())
    {
        rhoUf_() = fvc::interpolate(rho*U);
        surfaceVectorField n(mesh.Sf()/mesh.magSf());
        rhoUf_() += n*(fvc::absolute(phi, rho,U)/mesh.magSf() - (n & rhoUf_()));

    }
    