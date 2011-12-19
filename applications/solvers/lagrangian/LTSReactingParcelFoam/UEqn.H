    fvVectorMatrix UEqn
    (
        //pZones.ddt(rho, U)
        fvm::ddt(rho, U)
      + fvm::div(phi, U)
      + turbulence->divDevRhoReff(U)
     ==
        rho.dimensionedInternalField()*g
      + parcels.SU(U)
      + sources(rho, U)
    );

    sources.constrain(UEqn);

    pZones.addResistance(UEqn);

    if (pimple.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p));
    }
