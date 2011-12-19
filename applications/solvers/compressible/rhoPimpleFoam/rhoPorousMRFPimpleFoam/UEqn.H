// Solve the Momentum equation

tmp<fvVectorMatrix> UEqn
(
    //pZones.ddt(rho, U)
    fvm::ddt(rho, U)
  + fvm::div(phi, U)
  + turbulence->divDevRhoReff(U)
);

UEqn().relax();

mrfZones.addCoriolis(rho, UEqn());
pZones.addResistance(UEqn());

volScalarField rAU(1.0/UEqn().A());

if (pimple.momentumPredictor())
{
    solve(UEqn() == -fvc::grad(p));
}
