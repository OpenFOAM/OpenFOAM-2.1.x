// Solve the Momentum equation

tmp<fvVectorMatrix> UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(phi, U)
  + turbulence->divDevRhoReff(U)
);

UEqn().relax();

volScalarField rAU(1.0/UEqn().A());

if (pimple.momentumPredictor())
{
    solve(UEqn() == -fvc::grad(p));
    K = 0.5*magSqr(U);
}
