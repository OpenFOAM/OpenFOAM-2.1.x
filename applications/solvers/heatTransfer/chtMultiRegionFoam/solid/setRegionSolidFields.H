    fvMesh& mesh = solidRegions[i];
    basicSolidThermo& thermo = thermos[i];

    tmp<volScalarField> trho = thermo.rho();
    const volScalarField& rho = trho();

    tmp<volScalarField> tcp = thermo.Cp();
    const volScalarField& cp = tcp();

    tmp<volScalarField> tkappa = thermo.K();
    const volScalarField& kappa = tkappa();
    //tmp<volSymmTensorField> tkappa = thermo.directionalK();
    //const volSymmTensorField& kappa = tkappa();

    volScalarField& T = thermo.T();
