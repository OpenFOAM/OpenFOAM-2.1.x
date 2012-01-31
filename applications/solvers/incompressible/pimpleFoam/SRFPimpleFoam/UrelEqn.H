    // Relative momentum predictor
    tmp<fvVectorMatrix> UrelEqn
    (
        fvm::ddt(Urel)
      + fvm::div(phi, Urel)
      + turbulence->divDevReff(Urel)
      + SRF->Su()
    );

    UrelEqn().relax();

    sources.constrain(UrelEqn());

    solve(UrelEqn() == -fvc::grad(p) + sources(Urel));
