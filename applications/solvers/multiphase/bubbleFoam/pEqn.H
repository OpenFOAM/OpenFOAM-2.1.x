{
    surfaceScalarField alpha1f(fvc::interpolate(alpha1));
    surfaceScalarField alpha2f(scalar(1) - alpha1f);

    volScalarField rAU1(1.0/U1Eqn.A());
    volScalarField rAU2(1.0/U2Eqn.A());

    surfaceScalarField rAU1f(fvc::interpolate(rAU1));
    surfaceScalarField rAU2f(fvc::interpolate(rAU2));

    volVectorField HbyA1("HbyA1", U1);
    HbyA1 = rAU1*U1Eqn.H();

    volVectorField HbyA2("HbyA2", U2);
    HbyA2 = rAU2*U2Eqn.H();

    surfaceScalarField phiDrag1
    (
        fvc::interpolate(alpha2/rho1*dragCoef*rAU1)*phi2
      + rAU1f*(g & mesh.Sf())
    );
    surfaceScalarField phiDrag2
    (
        fvc::interpolate(alpha1/rho2*dragCoef*rAU2)*phi1
      + rAU2f*(g & mesh.Sf())
    );

    forAll(p.boundaryField(), patchi)
    {
        if (isA<zeroGradientFvPatchScalarField>(p.boundaryField()[patchi]))
        {
            phiDrag1.boundaryField()[patchi] = 0.0;
            phiDrag2.boundaryField()[patchi] = 0.0;
        }
    }

    surfaceScalarField phiHbyA1
    (
        (fvc::interpolate(HbyA1) & mesh.Sf())
      + fvc::ddtPhiCorr(rAU1, U1, phi1)
      + phiDrag1
    );

    surfaceScalarField phiHbyA2
    (
        (fvc::interpolate(HbyA2) & mesh.Sf())
      + fvc::ddtPhiCorr(rAU2, U2, phi2)
      + phiDrag2
    );

    surfaceScalarField phiHbyA("phiHbyA", alpha1f*phiHbyA1 + alpha2f*phiHbyA2);

    surfaceScalarField Dp
    (
        "Dp",
        alpha1f*rAU1f/rho1 + alpha2f*rAU2f/rho2
    );

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(Dp, p) == fvc::div(phiHbyA)
        );

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            surfaceScalarField SfGradp(pEqn.flux()/Dp);

            phi1 = phiHbyA1 - rAU1f*SfGradp/rho1;
            phi2 = phiHbyA2 - rAU2f*SfGradp/rho2;
            phi = alpha1f*phi1 + alpha2f*phi2;

            p.relax();
            SfGradp = pEqn.flux()/Dp;

            U1 = HbyA1 + (fvc::reconstruct(phiDrag1 - rAU1f*SfGradp/rho1));
            U1.correctBoundaryConditions();

            U2 = HbyA2 + (fvc::reconstruct(phiDrag2 - rAU2f*SfGradp/rho2));
            U2.correctBoundaryConditions();

            U = alpha1*U1 + alpha2*U2;
        }
    }
}

#include "continuityErrs.H"
