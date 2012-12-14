    // Construct the Momentum equation

    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff(U)
      ==
        sources(U)
    );

    mrfZones.addCoriolis(UEqn());

    UEqn().relax();

    sources.constrain(UEqn());

    // Include the porous media resistance and solve the momentum equation
    // either implicit in the tensorial resistance or transport using by
    // including the spherical part of the resistance in the momentum diagonal

    tmp<volScalarField> trAU;
    tmp<volTensorField> trTU;

    if (pressureImplicitPorosity)
    {
        tmp<volTensorField> tTU = tensor(I)*UEqn().A();
        pZones.addResistance(UEqn(), tTU());
        trTU = inv(tTU());
        trTU().rename("rAU");

        volVectorField gradp(fvc::grad(p));

        for (int UCorr=0; UCorr<nUCorr; UCorr++)
        {
            U = trTU() & (UEqn().H() - gradp);
        }
        U.correctBoundaryConditions();
    }
    else
    {
        pZones.addResistance(UEqn());

        solve(UEqn() == -fvc::grad(p));

        trAU = 1.0/UEqn().A();
        trAU().rename("rAU");
    }
