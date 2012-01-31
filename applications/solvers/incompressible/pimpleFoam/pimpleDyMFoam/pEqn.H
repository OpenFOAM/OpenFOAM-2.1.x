U = rAU*(UEqn() == sources(U))().H();

if (pimple.nCorrPISO() <= 1)
{
    UEqn.clear();
}

phi = (fvc::interpolate(U) & mesh.Sf());

if (ddtPhiCorr)
{
    phi += fvc::ddtPhiCorr(rAU, U, phi);
}

if (p.needReference())
{
    fvc::makeRelative(phi, U);
    adjustPhi(phi, U, p);
    fvc::makeAbsolute(phi, U);
}

while (pimple.correctNonOrthogonal())
{
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

// Make the fluxes relative to the mesh motion
fvc::makeRelative(phi, U);

U -= rAU*fvc::grad(p);
U.correctBoundaryConditions();
sources.correct(U);
