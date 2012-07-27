    volVectorField Ur(U1 - U2);
    volScalarField magUr(mag(Ur));

    volScalarField Cd1MagUr
    (
        (24.0*nu2/d1)*(scalar(1) + 0.15*pow(d1*magUr/nu2, 0.687))
    );

    volScalarField Cd2MagUr
    (
        (24.0*nu1/d2)*(scalar(1) + 0.15*pow(d2*magUr/nu1, 0.687))
    );

    volScalarField dragCoef
    (
        "Cd",
        0.75*(alpha2*rho2*Cd1MagUr/d1 + alpha1*rho1*Cd2MagUr/d2)
    );

    volVectorField liftCoeff
    (
        Cl*(alpha2*rho2 + alpha1*rho1)*(Ur ^ fvc::curl(U))
    );
