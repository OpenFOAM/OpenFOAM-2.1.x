{
    volScalarField rAU(1.0/UEqn.A());
    surfaceScalarField rAUf(fvc::interpolate(rAU));

    tmp<fvScalarMatrix> p_rghEqnComp;

    if (pimple.transonic())
    {
        p_rghEqnComp =
        (
            fvm::ddt(p_rgh)
          + fvm::div(phi, p_rgh)
          - fvm::Sp(fvc::div(phi), p_rgh)
        );
    }
    else
    {
        p_rghEqnComp =
        (
            fvm::ddt(p_rgh)
          + fvc::div(phi, p_rgh)
          - fvc::Sp(fvc::div(phi), p_rgh)
        );
    }


    U = rAU*UEqn.H();

    surfaceScalarField phiU
    (
        "phiU",
        (fvc::interpolate(U) & mesh.Sf())
      + fvc::ddtPhiCorr(rAU, rho, U, phi)
    );

    phi = phiU +
        (
            fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)
          - ghf*fvc::snGrad(rho)
        )*rAUf*mesh.magSf();

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqnIncomp
        (
            fvc::div(phi)
          - fvm::laplacian(rAUf, p_rgh)
        );

        solve
        (
            (
                max(alpha1, scalar(0))*(psi1/rho1)
              + max(alpha2, scalar(0))*(psi2/rho2)
            )
           *p_rghEqnComp()
          + p_rghEqnIncomp,
            mesh.solver(p_rgh.select(pimple.finalInnerIter()))
        );

        if (pimple.finalNonOrthogonalIter())
        {
            dgdt =
                (pos(alpha2)*(psi2/rho2) - pos(alpha1)*(psi1/rho1))
               *(p_rghEqnComp & p_rgh);
            phi += p_rghEqnIncomp.flux();
        }
    }

    U += rAU*fvc::reconstruct((phi - phiU)/rAUf);
    U.correctBoundaryConditions();

    p = max
    (
        (p_rgh + gh*(alpha1*rho10 + alpha2*rho20))
       /(1.0 - gh*(alpha1*psi1 + alpha2*psi2)),
        pMin
    );

    rho1 = rho10 + psi1*p;
    rho2 = rho20 + psi2*p;

    Info<< "max(U) " << max(mag(U)).value() << endl;
    Info<< "min(p_rgh) " << min(p_rgh).value() << endl;

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi, U);
}
