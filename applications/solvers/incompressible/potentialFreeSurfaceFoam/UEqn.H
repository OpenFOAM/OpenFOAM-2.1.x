tmp<fvVectorMatrix> UEqn
(
    fvm::ddt(U)
  + fvm::div(phi, U)
  + turbulence->divDevReff(U)
  ==
    sources(U)
);


UEqn().relax();

sources.constrain(UEqn());

if (pimple.momentumPredictor())
{
    solve(UEqn() == -fvc::grad(p_gh));
}
