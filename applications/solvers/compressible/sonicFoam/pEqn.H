rho = thermo.rho();

volScalarField rAU(1.0/UEqn.A());
U = rAU*UEqn.H();

surfaceScalarField phid
(
    "phid",
    fvc::interpolate(psi)
   *(
        (fvc::interpolate(U) & mesh.Sf())
      + fvc::ddtPhiCorr(rAU, rho, U, phi)
    )
);

for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
{
    fvScalarMatrix pEqn
    (
        fvm::ddt(psi, p)
      + fvm::div(phid, p)
      - fvm::laplacian(rho*rAU, p)
    );

    pEqn.solve();

    if (nonOrth == nNonOrthCorr)
    {
        phi = pEqn.flux();
    }
}

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

U -= rAU*fvc::grad(p);
U.correctBoundaryConditions();
