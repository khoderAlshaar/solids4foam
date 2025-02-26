{
    fvScalarMatrix EEqn
    (
        fvm::ddt(rho_, e_) + fvm::div(phi(), e_)
      + fvc::ddt(rho_, K_) + fvc::div(phi(), K_)
      + fvc::div(fvc::absolute(phi()/fvc::interpolate(rho_), U()), p(), "div(phiv,p)")
      - fvm::laplacian(turbulence_->alphaEff(), e_)
     ==
        options()(rho_, e_)
    );

    EEqn.relax();

    options().constrain(EEqn);

    EEqn.solve();

    options().correct(e_);

    thermo_.correct();
}
