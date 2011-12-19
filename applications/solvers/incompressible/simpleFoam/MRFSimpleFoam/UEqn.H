    // Momentum predictor

    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff(U)
      ==
        sources(U)
    );

    mrfZones.addCoriolis(UEqn());

    sources.constrain(UEqn());

    UEqn().relax();

    solve(UEqn() == -fvc::grad(p));
