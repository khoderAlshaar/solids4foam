// Solve the Momentum equation

MRF.correctBoundaryVelocity(U());

fvVectorMatrix UEqn
(
    fvm::ddt(rho_, U()) + fvm::div(phi(), U())
  + MRF.DDt(rho_, U())
  + turbulence_->divDevRhoReff(U())
 ==
    options()(rho_, U())
);

UEqn.relax();

options().constrain(UEqn);

if (pimple().momentumPredictor())
{
    solve(UEqn == -fvc::grad(p()));

    options().correct(U());
    K_ = 0.5*magSqr(U());
}
