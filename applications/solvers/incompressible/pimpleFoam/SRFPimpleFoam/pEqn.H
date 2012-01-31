volScalarField rAUrel(1.0/UrelEqn().A());
Urel = rAUrel*(UrelEqn() == sources(Urel))().H();

if (pimple.nCorrPISO() <= 1)
{
    UrelEqn.clear();
}

phi = (fvc::interpolate(Urel) & mesh.Sf())
    + fvc::ddtPhiCorr(rAUrel, Urel, phi);

adjustPhi(phi, Urel, p);

// Non-orthogonal pressure corrector loop
while (pimple.correctNonOrthogonal())
{
    // Pressure corrector
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAUrel, p) == fvc::div(phi)
    );

    pEqn.setReference(pRefCell, pRefValue);

    pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

    if (pimple.finalNonOrthogonalIter())
    {
        phi -= pEqn.flux();
    }
}

#include "continuityErrs.H"

p.relax();

// Momentum corrector
Urel -= rAUrel*fvc::grad(p);
Urel.correctBoundaryConditions();
sources.correct(Urel);
