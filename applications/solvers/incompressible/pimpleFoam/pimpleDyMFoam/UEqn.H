// Solve the Momentum equation

tmp<fvVectorMatrix> UEqn
(
    fvm::ddt(U)
  + fvm::div(phi, U)
  + turbulence->divDevReff(U)
);

UEqn().relax();

rAU = 1.0/UEqn().A();

if (pimple.momentumPredictor())
{
    solve(UEqn() == -fvc::grad(p));
}
