{
    rho1 = rho10 + psi1*p;
    rho2 = rho20 + psi2*p;

    tmp<fvScalarMatrix> pEqnComp1;
    tmp<fvScalarMatrix> pEqnComp2;

    //if (transonic)
    //{
    //}
    //else
    {
        surfaceScalarField phid1("phid1", fvc::interpolate(psi1)*phi1);
        surfaceScalarField phid2("phid2", fvc::interpolate(psi2)*phi2);

        pEqnComp1 =
            fvc::ddt(rho1) + psi1*correction(fvm::ddt(p))
          + fvc::div(phid1, p)
          - fvc::Sp(fvc::div(phid1), p);

        pEqnComp2 =
            fvc::ddt(rho2) + psi2*correction(fvm::ddt(p))
          + fvc::div(phid2, p)
          - fvc::Sp(fvc::div(phid2), p);
    }

    surfaceScalarField alpha1f(fvc::interpolate(alpha1));
    surfaceScalarField alpha2f(scalar(1) - alpha1f);

    volVectorField U10 = U1;
    volVectorField U20 = U2;

    rAU1 = 1.0/U1Eqn.A();
    volScalarField rAU2(1.0/U2Eqn.A());

    surfaceScalarField rAlphaAU1f = fvc::interpolate(alpha1*rAU1);
    surfaceScalarField rAlphaAU2f = fvc::interpolate(alpha2*rAU2);

    U1 = rAU1*U1Eqn.H();
    U2 = rAU2*U2Eqn.H();

    {
        surfaceScalarField phi10 = phi1;

        phi1 =
        (
            (fvc::interpolate(U1) & mesh.Sf())
          + fvc::ddtPhiCorr(rAU1, alpha1, U1, phi1)
          + fvc::interpolate((1.0/rho1)*rAU1*dragCoeff)*phi2
          + rAlphaAU1f*(g & mesh.Sf())
        );

        phi2 =
        (
            (fvc::interpolate(U2) & mesh.Sf())
          + fvc::ddtPhiCorr(rAU2, alpha2, U2, phi2)
          + fvc::interpolate((1.0/rho2)*rAU2*dragCoeff)*phi10
          + rAlphaAU2f*(g & mesh.Sf())
        );
    }

    if (g0.value() > 0.0)
    {
        phi1 -= ppMagf*fvc::snGrad(alpha1)*mesh.magSf();
    }

    if (kineticTheory.on())
    {
        phi1 -=
            fvc::interpolate((1.0/rho1)*rAU1)
           *fvc::snGrad(kineticTheory.pa())*mesh.magSf();
    }

    phi = alpha1f*phi1 + alpha2f*phi2;
    adjustPhi(phi, U, p);

    surfaceScalarField Dp
    (
        "(rho*1|A(U))",
        mag(alpha1f*rAlphaAU1f/fvc::interpolate(rho1)
      + alpha2f*rAlphaAU2f/fvc::interpolate(rho2))
    );

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqnIncomp
        (
            fvc::div(phi)
          - fvm::laplacian(Dp, p)
        );

        solve
        (
            (
                (alpha1/rho1)*pEqnComp1()
              + (alpha2/rho2)*pEqnComp2()
            )
          + pEqnIncomp,
            mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            surfaceScalarField mSfGradp = pEqnIncomp.flux()/Dp;
            phi1 += rAlphaAU1f*mSfGradp/fvc::interpolate(rho1);
            phi2 += rAlphaAU2f*mSfGradp/fvc::interpolate(rho2);
            phi = alpha1f*phi1 + alpha2f*phi2;

            dgdt =
            (
                pos(alpha2)*(pEqnComp2 & p)/rho2
              - pos(alpha1)*(pEqnComp1 & p)/rho1
            );

            p.relax();
            mSfGradp = pEqnIncomp.flux()/Dp;

            U1 += (1.0/rho1)*rAU1*dragCoeff*U20
               + fvc::reconstruct
                 (
                     rAlphaAU1f*(g & mesh.Sf())
                   + rAlphaAU1f*mSfGradp/fvc::interpolate(rho1)
                 );
            U1.correctBoundaryConditions();

            U2 += (1.0/rho2)*rAU2*dragCoeff*U10
               + fvc::reconstruct
                 (
                     rAlphaAU2f*(g & mesh.Sf())
                   + rAlphaAU2f*mSfGradp/fvc::interpolate(rho2)
                 );
            U2.correctBoundaryConditions();

            U = alpha1*U1 + alpha2*U2;
        }
    }

    p = max(p, pMin);

    rho1 = rho10 + psi1*p;
    rho2 = rho20 + psi2*p;

    K1 = 0.5*magSqr(U1);
    K2 = 0.5*magSqr(U1);

    dpdt = fvc::ddt(p);
}
