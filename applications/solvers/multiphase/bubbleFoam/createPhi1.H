    IOobject phi1Header
    (
        "phi1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );

    autoPtr<surfaceScalarField> phi1Ptr(NULL);

    if (phi1Header.headerOk())
    {
        Info<< "Reading face flux field phi1\n" << endl;

        phi1Ptr.reset
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "phi1",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );
    }
    else
    {
        Info<< "Calculating face flux field phi1\n" << endl;

        wordList phiTypes
        (
            U1.boundaryField().size(),
            calculatedFvPatchScalarField::typeName
        );

        forAll(U1.boundaryField(), i)
        {
            if (isA<fixedValueFvPatchVectorField>(U1.boundaryField()[i]))
            {
                phiTypes[i] = fixedValueFvPatchScalarField::typeName;
            }
        }

        phi1Ptr.reset
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "phi1",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fvc::interpolate(U1) & mesh.Sf(),
                phiTypes
            )
        );
    }

    surfaceScalarField& phi1 = phi1Ptr();
