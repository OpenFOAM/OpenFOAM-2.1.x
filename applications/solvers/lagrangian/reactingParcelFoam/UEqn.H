    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(phi, U)
      + turbulence->divDevRhoReff(U)
     ==
        rho.dimensionedInternalField()*g
      + parcels.SU(U)
    );

    UEqn.relax();

    if (pimple.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p) + sources(rho, U));
        K = 0.5*magSqr(U);
    }
