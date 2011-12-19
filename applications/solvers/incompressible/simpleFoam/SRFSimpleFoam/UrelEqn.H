    // Relative momentum predictor

    tmp<fvVectorMatrix> UrelEqn
    (
        fvm::div(phi, Urel)
      + turbulence->divDevReff(Urel)
      + SRF->Su()
     ==
        sources(Urel)
    );

    UrelEqn().relax();

    sources.constrain(UrelEqn());

    solve(UrelEqn() == -fvc::grad(p));
