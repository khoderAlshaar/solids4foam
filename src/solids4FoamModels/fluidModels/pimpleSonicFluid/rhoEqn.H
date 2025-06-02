    // fvScalarMatrix rhoEqn
    // (
    //     fvm::ddt(rho)
    //   + fvc::div(phi)
    //  ==
    //     options()(rho)
    // );

    // options().constrain(rhoEqn);

    // rhoEqn.solve();

    // options().correct(rho);
{
    solve(fvm::ddt(rho) + fvc::div(phi));
}
