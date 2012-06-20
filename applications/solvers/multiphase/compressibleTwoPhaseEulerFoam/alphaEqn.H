surfaceScalarField alphaPhi1("alphaPhi1", phi1);
surfaceScalarField alphaPhi2("alphaPhi2", phi2);

{
    word scheme("div(phi,alpha1)");
    word schemer("div(phir,alpha1)");

    surfaceScalarField phic("phic", phi);
    surfaceScalarField phir("phir", phi1 - phi2);

    if (g0.value() > 0.0)
    {
        surfaceScalarField alpha1f(fvc::interpolate(alpha1));
        surfaceScalarField phipp(ppMagf*fvc::snGrad(alpha1)*mesh.magSf());
        phir += phipp;
        phic += fvc::interpolate(alpha1)*phipp;
    }

    for (int acorr=0; acorr<nAlphaCorr; acorr++)
    {
        volScalarField::DimensionedInternalField Sp
        (
            IOobject
            (
                "Sp",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("Sp", dgdt.dimensions(), 0.0)
        );

        volScalarField::DimensionedInternalField Su
        (
            IOobject
            (
                "Su",
                runTime.timeName(),
                mesh
            ),
            // Divergence term is handled explicitly to be
            // consistent with the explicit transport solution
            fvc::div(phi)*min(alpha1, scalar(1))
        );

        forAll(dgdt, celli)
        {
            if (dgdt[celli] > 0.0 && alpha1[celli] > 0.0)
            {
                Sp[celli] -= dgdt[celli]*alpha1[celli];
                Su[celli] += dgdt[celli]*alpha1[celli];
            }
            else if (dgdt[celli] < 0.0 && alpha1[celli] < 1.0)
            {
                Sp[celli] += dgdt[celli]*(1.0 - alpha1[celli]);
            }
        }


        fvScalarMatrix alpha1Eqn
        (
             fvm::ddt(alpha1)
           + fvm::div(phic, alpha1, scheme)
           + fvm::div(-fvc::flux(-phir, alpha2, schemer), alpha1, schemer)
          ==
             fvm::Sp(Sp, alpha1) + Su
        );

        if (g0.value() > 0.0)
        {
            ppMagf =
                fvc::interpolate((1.0/rho1)*rAU1)
               *fvc::interpolate
                (
                    g0*min(exp(preAlphaExp*(alpha1 - alphaMax)), expMax)
                );

            alpha1Eqn -= fvm::laplacian
            (
                (fvc::interpolate(alpha1) + scalar(0.0001))*ppMagf,
                alpha1,
                "laplacian(alphaPpMag,alpha1)"
            );
        }

        alpha1Eqn.relax();
        alpha1Eqn.solve();

        //***HGW temporary boundedness-fix pending the introduction of MULES
        alpha1 = max(min(alpha1, scalar(1)), scalar(0));

        #include "packingLimiter.H"

        alphaPhi1 = alpha1Eqn.flux();
        alphaPhi2 = phi - alphaPhi1;
        alpha2 = scalar(1) - alpha1;

        Info<< "Dispersed phase volume fraction = "
            << alpha1.weightedAverage(mesh.V()).value()
            << "  Min(alpha1) = " << min(alpha1).value()
            << "  Max(alpha1) = " << max(alpha1).value()
            << endl;
    }
}

rho = alpha1*rho1 + alpha2*rho2;
