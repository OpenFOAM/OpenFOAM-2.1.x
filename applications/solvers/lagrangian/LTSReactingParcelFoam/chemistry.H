if (chemistry.chemistry())
{
    Info<< "Solving chemistry" << endl;

    // update reaction rates
    chemistry.calculate();

    // turbulent time scale
    if (turbulentReaction)
    {
        typedef DimensionedField<scalar, volMesh> dsfType;

        const dimensionedScalar e0("e0", sqr(dimLength)/pow3(dimTime), SMALL);

        const dsfType tk =
            Cmix*sqrt(turbulence->muEff()/rho/(turbulence->epsilon() + e0));

        const dsfType tc = chemistry.tc()().dimensionedInternalField();

        kappa = tc/(tc + tk);
    }
    else
    {
        kappa = 1.0;
    }

    chemistrySh = kappa*chemistry.Sh()();
}
