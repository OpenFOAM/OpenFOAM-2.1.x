volScalarField rAU(1.0/UEqn.A());

surfaceScalarField rAUf
(
    "(rho*(1|A(U)))",
    fvc::interpolate(rho)*fvc::interpolate(rAU)
);

U = rAU*UEqn.H();
phi =
    fvc::interpolate(rho)
   *(
       (fvc::interpolate(U) & mesh.Sf())
     + fvc::ddtPhiCorr(rAU, rho, U, phi)
    );

surfaceScalarField phiU("phiU", phi);
phi -= ghf*fvc::snGrad(rho)*rAUf*mesh.magSf();

while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix p_rghEqn
    (
        fvm::laplacian(rAUf, p_rgh) == fvc::ddt(rho) + fvc::div(phi)
    );

    p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

    p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

    if (pimple.finalNonOrthogonalIter())
    {
        phi -= p_rghEqn.flux();
    }
}

p == p_rgh + rho*gh;

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

U += rAU*fvc::reconstruct((phi - phiU)/rAUf);
U.correctBoundaryConditions();
