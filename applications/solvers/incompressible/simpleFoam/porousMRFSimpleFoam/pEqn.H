if (pressureImplicitPorosity)
{
    U = trTU()&UEqn().H();
}
else
{
    U = trAU()*UEqn().H();
}

UEqn.clear();
phi = fvc::interpolate(U) & mesh.Sf();
mrfZones.relativeFlux(phi);
adjustPhi(phi, U, p);

while (simple.correctNonOrthogonal())
{
    tmp<fvScalarMatrix> tpEqn;

    if (pressureImplicitPorosity)
    {
        tpEqn = (fvm::laplacian(trTU(), p) == fvc::div(phi));
    }
    else
    {
        tpEqn = (fvm::laplacian(trAU(), p) == fvc::div(phi));
    }

    tpEqn().setReference(pRefCell, pRefValue);

    tpEqn().solve();

    if (simple.finalNonOrthogonalIter())
    {
        phi -= tpEqn().flux();
    }
}

#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

if (pressureImplicitPorosity)
{
    U -= trTU()&fvc::grad(p);
}
else
{
    U -= trAU()*fvc::grad(p);
}

U.correctBoundaryConditions();
