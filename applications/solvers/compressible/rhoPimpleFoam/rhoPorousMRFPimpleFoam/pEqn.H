rho = thermo.rho();
rho = max(rho, rhoMin);
rho = min(rho, rhoMax);
rho.relax();

volScalarField rAU(1.0/UEqn().A());
U = rAU*UEqn().H();

if (pimple.nCorrPISO() <= 1)
{
    UEqn.clear();
}

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
    mrfZones.relativeFlux(fvc::interpolate(psi), phid);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi, p)
          + fvm::div(phid, p)
          - fvm::laplacian(rho*rAU, p)
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
    mrfZones.relativeFlux(fvc::interpolate(rho), phi);

    while (pimple.correctNonOrthogonal())
    {
        // Pressure corrector
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi, p)
          + fvc::div(phi)
          - fvm::laplacian(rho*rAU, p)
        );

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi += pEqn.flux();
        }
    }
}

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

rho = thermo.rho();
rho = max(rho, rhoMin);
rho = min(rho, rhoMax);
rho.relax();
Info<< "rho max/min : " << max(rho).value()
    << " " << min(rho).value() << endl;

U -= rAU*fvc::grad(p);
U.correctBoundaryConditions();
K = 0.5*magSqr(U);

dpdt = fvc::ddt(p);
