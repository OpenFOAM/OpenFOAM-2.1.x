    // Momentum predictor

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

    solve(UEqn() == -fvc::grad(p));
