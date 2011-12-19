    Info<< "Creating combustion model\n" << endl;

    autoPtr<combustionModels::rhoChemistryCombustionModel> combustion
    (
        combustionModels::rhoChemistryCombustionModel::New
        (
            mesh
        )
    );

    rhoChemistryModel& chemistry = combustion->pChemistry();

    hsReactionThermo& thermo = chemistry.thermo();

    SLGThermo slgThermo(mesh, thermo);

    basicMultiComponentMixture& composition = thermo.composition();
    PtrList<volScalarField>& Y = composition.Y();

    const word inertSpecie(thermo.lookup("inertSpecie"));

    if (!composition.contains(inertSpecie))
    {
        FatalErrorIn(args.executable())
            << "Specified inert specie '" << inertSpecie << "' not found in "
            << "species list. Available species:" << composition.species()
            << exit(FatalError);
    }

    volScalarField& p = thermo.p();
    volScalarField& hs = thermo.hs();
    const volScalarField& T = thermo.T();
    const volScalarField& psi = thermo.psi();

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo.rho()
    );

    Info<< "\nReading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    #include "compressibleCreatePhi.H"

    dimensionedScalar rhoMax
    (
        mesh.solutionDict().subDict("PIMPLE").lookup("rhoMax")
    );

    dimensionedScalar rhoMin
    (
        mesh.solutionDict().subDict("PIMPLE").lookup("rhoMin")
    );

    Info<< "Creating turbulence model\n" << endl;
    autoPtr<compressible::turbulenceModel> turbulence
    (
        compressible::turbulenceModel::New
        (
            rho,
            U,
            phi,
            thermo
        )
    );

    // Set the turbulence into the combustion model
    combustion->setTurbulence(turbulence());

    Info<< "Creating multi-variate interpolation scheme\n" << endl;
    multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

    forAll(Y, i)
    {
        fields.add(Y[i]);
    }
    fields.add(hs);

    volScalarField dQ
    (
        IOobject
        (
            "dQ",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("dQ", dimEnergy/dimTime, 0.0)
    );


    volScalarField rDeltaT
    (
        IOobject
        (
            "rDeltaT",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("one", dimless/dimTime, 1),
        zeroGradientFvPatchScalarField::typeName
    );
