{
    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, e) + fvm::div(phi, e)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + fvc::div(fvc::absolute(phi/fvc::interpolate(rho), U), p, "div(phiv,p)")
      - fvm::laplacian(turbulence_->alphaEff(), e)
     ==
        options()(rho, e)
    );

    EEqn.relax();

    options().constrain(EEqn);

    EEqn.solve();

    options().correct(e);

    thermo_.correct();
}

