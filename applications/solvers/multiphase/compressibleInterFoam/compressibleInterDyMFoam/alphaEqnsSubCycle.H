{
    label nAlphaCorr(readLabel(pimple.dict().lookup("nAlphaCorr")));

    label nAlphaSubCycles(readLabel(pimple.dict().lookup("nAlphaSubCycles")));

    surfaceScalarField phic(mag(phi/mesh.magSf()));
    phic = min(interface.cAlpha()*phic, max(phic));

    fvc::makeAbsolute(phi, U);
    volScalarField divU(fvc::div(phi));
    fvc::makeRelative(phi, U);

    if (nAlphaSubCycles > 1)
    {
        dimensionedScalar totalDeltaT = runTime.deltaT();
        surfaceScalarField rhoPhiSum(0.0*rhoPhi);

        for
        (
            subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
            !(++alphaSubCycle).end();
        )
        {
            #include "alphaEqns.H"
            rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
        }

        rhoPhi = rhoPhiSum;
    }
    else
    {
        #include "alphaEqns.H"
    }

    // correct interface on first PIMPLE corrector
    if (pimple.corr() == 1)
    {
        interface.correct();
    }
}
