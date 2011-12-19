if (turbulence)
{
    if (mesh.changing())
    {
        y.correct();
    }

    dimensionedScalar kMin("kMin", k.dimensions(), SMALL);
    dimensionedScalar epsilonMin("epsilonMin", epsilon.dimensions(), SMALL);

    volScalarField divU(fvc::div(phi/fvc::interpolate(rho)));

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField G(2*mut*(tgradU() && dev(symm(tgradU()))));
    tgradU.clear();

    volScalarField Gcoef
    (
        Cmu*k/sigmak*(g & fvc::grad(rho))/(epsilon + epsilonMin)
    );

    #include "wallFunctions.H"

    // Dissipation equation
    fvScalarMatrix epsEqn
    (
        fvm::ddt(rho, epsilon)
      + fvm::div(phi, epsilon)
      - fvm::laplacian
        (
            mut/sigmaEps + mul, epsilon,
            "laplacian(DepsilonEff,epsilon)"
        )
     ==
        C1*G*epsilon/k
      - fvm::SuSp(C1*(1.0 - C3)*Gcoef + (2.0/3.0*C1)*rho*divU, epsilon)
      - fvm::Sp(C2*rho*epsilon/k, epsilon)
    );

    #include "wallDissipation.H"

    epsEqn.relax();
    epsEqn.solve();

    bound(epsilon, epsilonMin);


    // Turbulent kinetic energy equation
    fvScalarMatrix kEqn
    (
        fvm::ddt(rho, k)
      + fvm::div(phi, k)
      - fvm::laplacian
        (
            mut/sigmak + mul, k,
            "laplacian(DkEff,k)"
        )
     ==
        G
      - fvm::SuSp(Gcoef + 2.0/3.0*rho*divU, k)
      - fvm::Sp(rho*epsilon/k, k)
    );

    kEqn.relax();
    kEqn.solve();

    bound(k, kMin);


    //- Re-calculate viscosity
    mut = rho*Cmu*sqr(k)/(epsilon + epsilonMin);

    #include "wallViscosity.H"
}

mu = mut + mul;
