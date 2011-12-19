Info<< "Reading thermophysicalProperties\n" << endl;

// Pr defined as a separate constant to enable calculation of k, currently
// inaccessible through thermo
IOdictionary thermophysicalProperties
(
    IOobject
    (
        "thermophysicalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar Pr
(
    dimensionedScalar::lookupOrDefault
    (
        "Pr",
        thermophysicalProperties,
        1.0
    )
);

