rho = thermo.rho();

volScalarField rAU(1.0/UEqn.A());
U = UEqn.H()/UEqn.A();

surfaceScalarField phid
(
    "phid",
    fvc::interpolate(psi)
   *(
        (fvc::interpolate(U) & mesh.Sf())
      - fvc::meshPhi(rho, U)
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

    phi = pEqn.flux();
}

#include "compressibleContinuityErrs.H"

U -= rAU*fvc::grad(p);
U.correctBoundaryConditions();
