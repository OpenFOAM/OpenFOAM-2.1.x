    fvVectorMatrix UEqn
    (
        //pZones.ddt(rho, U)
        fvm::ddt(rho, U)
      + fvm::div(phi, U)
      + turbulence->divDevRhoReff(U)
     ==
        rho.dimensionedInternalField()*g
      + parcels.SU(U)
    );

    UEqn.relax();

    sources.constrain(UEqn);

    pZones.addResistance(UEqn);

    if (pimple.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p) + sources(rho, U));
        K = 0.5*magSqr(U);
    }

