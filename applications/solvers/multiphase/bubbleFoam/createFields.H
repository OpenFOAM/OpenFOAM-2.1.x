    Info<< "Reading field alpha1\n" << endl;
    volScalarField alpha1
    (
        IOobject
        (
            "alpha1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField alpha2
    (
        IOobject
        (
            "alpha2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        scalar(1) - alpha1
        //,alpha1.boundaryField().types()
    );

    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field U1\n" << endl;
    volVectorField U1
    (
        IOobject
        (
            "U1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field U2\n" << endl;
    volVectorField U2
    (
        IOobject
        (
            "U2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1*U1 + alpha2*U2
    );


    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar rho1
    (
        transportProperties.lookup("rho1")
    );

    dimensionedScalar rho2
    (
        transportProperties.lookup("rho2")
    );

    dimensionedScalar nu1
    (
        transportProperties.lookup("nu1")
    );

    dimensionedScalar nu2
    (
        transportProperties.lookup("nu2")
    );

    dimensionedScalar d1
    (
        transportProperties.lookup("d1")
    );

    dimensionedScalar d2
    (
        transportProperties.lookup("d2")
    );

    dimensionedScalar Cvm
    (
        transportProperties.lookup("Cvm")
    );

    dimensionedScalar Cl
    (
        transportProperties.lookup("Cl")
    );

    dimensionedScalar Ct
    (
        transportProperties.lookup("Ct")
    );

    #include "createPhi1.H"
    #include "createPhi2.H"

    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh
        ),
        fvc::interpolate(alpha1)*phi1
      + fvc::interpolate(alpha2)*phi2
    );

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh
        ),
        alpha1*rho1 + alpha2*rho2
    );

    #include "createRASTurbulence.H"

    Info<< "Calculating field DDtU1 and DDtU2\n" << endl;

    volVectorField DDtU1
    (
        fvc::ddt(U1)
      + fvc::div(phi1, U1)
      - fvc::div(phi1)*U1
    );

    volVectorField DDtU2
    (
        fvc::ddt(U2)
      + fvc::div(phi2, U2)
      - fvc::div(phi2)*U2
    );


    Info<< "Calculating field g.h\n" << endl;
    volScalarField gh("gh", g & mesh.C());


    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("PIMPLE"), pRefCell, pRefValue);
