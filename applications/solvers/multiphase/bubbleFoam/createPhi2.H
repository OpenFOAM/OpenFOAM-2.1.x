    IOobject phi2Header
    (
        "phi2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );

    autoPtr<surfaceScalarField> phi2Ptr(NULL);

    if (phi2Header.headerOk())
    {
        Info<< "Reading face flux field phi2\n" << endl;

        phi2Ptr.reset
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "phi2",
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
        Info<< "Calculating face flux field phi2\n" << endl;

        wordList phiTypes
        (
            U2.boundaryField().size(),
            calculatedFvPatchScalarField::typeName
        );

        forAll(U2.boundaryField(), i)
        {
            if (isA<fixedValueFvPatchVectorField>(U2.boundaryField()[i]))
            {
                phiTypes[i] = fixedValueFvPatchScalarField::typeName;
            }
        }

        phi2Ptr.reset
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "phi2",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fvc::interpolate(U2) & mesh.Sf(),
                phiTypes
            )
        );
    }

    surfaceScalarField& phi2 = phi2Ptr();
