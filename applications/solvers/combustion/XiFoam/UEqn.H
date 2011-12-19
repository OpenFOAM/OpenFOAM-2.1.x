    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(phi, U)
      + turbulence->divDevRhoReff(U)
     ==
        rho*g
    );

    UEqn.relax();

    if (pimple.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p));
        K = 0.5*magSqr(U);
    }
