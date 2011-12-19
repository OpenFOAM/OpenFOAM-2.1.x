volScalarField dragCoeff
(
    IOobject
    (
        "dragCoeff",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("dragCoeff", dimensionSet(1, -3, -1, 0, 0), 0)
);

volVectorField liftForce
(
    IOobject
    (
        "liftForce",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedVector("liftForce", dimensionSet(1, -2, -2, 0, 0), vector::zero)
);

volScalarField heatTransferCoeff
(
    IOobject
    (
        "heatTransferCoeff",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("heatTransferCoeff", dimensionSet(1, -1, -3, -1, 0), 0)
);

{
    volVectorField Ur = U1 - U2;
    volScalarField magUr = mag(Ur);

    if (dispersedPhase == "1")
    {
        dragCoeff = drag1->K(magUr);
        heatTransferCoeff = heatTransfer1->K(magUr);
    }
    else if (dispersedPhase == "2")
    {
        dragCoeff = drag2->K(magUr);
        heatTransferCoeff = heatTransfer2->K(magUr);
    }
    else if (dispersedPhase == "both")
    {
        dragCoeff =
        (
            alpha2*drag1->K(magUr)
          + alpha1*drag2->K(magUr)
        );

        heatTransferCoeff =
        (
            alpha2*heatTransfer1->K(magUr)
          + alpha1*heatTransfer2->K(magUr)
        );
    }
    else
    {
        FatalErrorIn(args.executable())
            << "dispersedPhase: " << dispersedPhase << " is incorrect"
            << exit(FatalError);
    }

    volScalarField alphaCoeff
    (
        (alpha1 + minInterfaceAlpha)*(alpha2 + minInterfaceAlpha)
    );
    dragCoeff *= alphaCoeff;
    heatTransferCoeff *= alphaCoeff;

    liftForce = Cl*(alpha1*rho1 + alpha2*rho2)*(Ur ^ fvc::curl(U));
}
