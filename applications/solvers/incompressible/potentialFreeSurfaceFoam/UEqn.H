tmp<fvVectorMatrix> UEqn
(
    fvm::ddt(U)
  + fvm::div(phi, U)
  + turbulence->divDevReff(U)
);


UEqn().relax();

sources.constrain(UEqn());

if (pimple.momentumPredictor())
{
    solve(UEqn() == -fvc::grad(p_gh) + sources(U));
}
