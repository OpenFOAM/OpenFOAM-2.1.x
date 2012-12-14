volScalarField rAU(1.0/UEqn().A());
surfaceScalarField rAUf(rAU.name() + 'f', fvc::interpolate(rAU));

U = rAU*(UEqn() == sources(U))().H();

if (pimple.nCorrPISO() <= 1)
{
    UEqn.clear();
}

phi = (fvc::interpolate(U) & mesh.Sf())
    + fvc::ddtPhiCorr(rAU, U, phi);

adjustPhi(phi, U, p_gh);


// Non-orthogonal pressure corrector loop
while (pimple.correctNonOrthogonal())
{ 
    fvScalarMatrix p_ghEqn
    (
        fvm::laplacian(rAUf, p_gh) == fvc::div(phi)
    );

    p_ghEqn.setReference(p_ghRefCell, p_ghRefValue);

    p_ghEqn.solve(mesh.solver(p_gh.select(pimple.finalInnerIter())));

    if (pimple.finalNonOrthogonalIter())
    {
        phi -= p_ghEqn.flux();
    }
}

#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p_gh.relax();

p = p_gh + (g & mesh.C());

U -= rAU*fvc::grad(p_gh);
U.correctBoundaryConditions();
sources.correct(U);
