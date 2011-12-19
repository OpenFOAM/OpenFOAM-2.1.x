    // Momentum predictor

    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff(U)
      ==
        sources(U)
    );

    UEqn().relax();

    sources.constrain(UEqn());

    solve(UEqn() == -fvc::grad(p));
