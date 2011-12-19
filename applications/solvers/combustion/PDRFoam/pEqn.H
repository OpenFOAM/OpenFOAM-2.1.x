rho = thermo.rho();

volScalarField rAU(1.0/UEqn.A());
U = invA & UEqn.H();

if (pimple.transonic())
{
    surfaceScalarField phid
    (
        "phid",
        fvc::interpolate(psi)
       *(
            (fvc::interpolate(U) & mesh.Sf())
          + fvc::ddtPhiCorr(rAU, rho, U, phi)
        )
    );

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            betav*fvm::ddt(psi, p)
          + fvm::div(phid, p)
          - fvm::laplacian(rho*invA, p)
        );

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi == pEqn.flux();
        }
    }
}
else
{
    phi =
        fvc::interpolate(rho)*
        (
            (fvc::interpolate(U) & mesh.Sf())
          + fvc::ddtPhiCorr(rAU, rho, U, phi)
        );

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            betav*fvm::ddt(psi, p)
          + fvc::div(phi)
          - fvm::laplacian(rho*invA, p)
        );

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi += pEqn.flux();
        }
    }
}

#include "rhoEqn.H"
#include "continuityErrs.H"

U -= invA & (betav*fvc::grad(p));
U.correctBoundaryConditions();
K = 0.5*magSqr(U);

dpdt = fvc::ddt(p);
