rho = thermo.rho();

volScalarField rAU(1.0/UEqn.A());
U = rAU*(UEqn == sources(rho, U))().H();

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
            fvm::ddt(psi, p)
          + fvm::div(phid, p)
          - fvm::laplacian(rho*rAU, p)
         ==
            coalParcels.Srho()
          + sources(psi, p, rho.name())
        );

        sources.constrain(pEqn, rho.name());

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
        fvc::interpolate(rho)
       *(
            (fvc::interpolate(U) & mesh.Sf())
          + fvc::ddtPhiCorr(rAU, rho, U, phi)
        );

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi, p)
          + fvc::div(phi)
          - fvm::laplacian(rho*rAU, p)
         ==
            coalParcels.Srho()
          + sources(psi, p, rho.name())
        );

        sources.constrain(pEqn, rho.name());

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi += pEqn.flux();
        }
    }
}

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

U -= rAU*fvc::grad(p);
U.correctBoundaryConditions();
sources.correct(U);

K = 0.5*magSqr(U);

dpdt = fvc::ddt(p);
