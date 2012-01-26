U = rAU*(UEqn() == sources(U))().H();

if (pimple.nCorrPISO() <= 1)
{
    UEqn.clear();
}

phi = (fvc::interpolate(U) & mesh.Sf())
    + fvc::ddtPhiCorr(rAU, U, phi);

adjustPhi(phi, U, p);

// Non-orthogonal pressure corrector loop
while (pimple.correctNonOrthogonal())
{
    // Pressure corrector
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAU, p) == fvc::div(phi)
    );

    pEqn.setReference(pRefCell, pRefValue);

    pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

    if (pimple.finalNonOrthogonalIter())
    {
        phi -= pEqn.flux();
    }
}

#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

U -= rAU*fvc::grad(p);
U.correctBoundaryConditions();
sources.correct(U);
