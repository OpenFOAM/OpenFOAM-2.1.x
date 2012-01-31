    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(phi, U)
      + turbulence->divDevRhoReff(U)
     ==
        rho.dimensionedInternalField()*g
      + coalParcels.SU(U)
      + limestoneParcels.SU(U)
    );

    UEqn.relax();

    sources.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p) + sources(rho, U));
        K = 0.5*magSqr(U);
    }
