{
    word scheme("div(phi,alpha1)");

    surfaceScalarField phir(phi1 - phi2);

    Info<< "Max Ur Courant Number = "
        << (
               max
               (
                   mesh.surfaceInterpolation::deltaCoeffs()*mag(phir)
                  /mesh.magSf()
                )*runTime.deltaT()
            ).value()
        << endl;

    for (int acorr=0; acorr<nAlphaCorr; acorr++)
    {
        fvScalarMatrix alpha1Eqn
        (
             fvm::ddt(alpha1)
           + fvm::div(phi, alpha1, scheme)
           + fvm::div(-fvc::flux(-phir, alpha2, scheme), alpha1, scheme)
        );
        alpha1Eqn.relax();
        alpha1Eqn.solve();

        /*
        fvScalarMatrix alpha2Eqn
        (
            fvm::ddt(alpha2)
          + fvm::div(phi, alpha2, scheme)
          + fvm::div
            (
                -fvc::flux(phir, scalar(1) - alpha2, scheme),
                alpha2,
                scheme
            )
        );
        alpha2Eqn.relax();
        alpha2Eqn.solve();

        alpha1 =
            0.5
           *(
                 scalar(1)
               + sqr(scalar(1) - alpha2)
               - sqr(scalar(1) - alpha1)
            );
        */

        alpha2 = scalar(1) - alpha1;
    }

    Info<< "Dispersed phase volume fraction = "
        << alpha1.weightedAverage(mesh.V()).value()
        << "  Min(alpha1) = " << min(alpha1).value()
        << "  Max(alpha1) = " << max(alpha1).value()
        << endl;
}

rho = alpha1*rho1 + alpha2*rho2;
