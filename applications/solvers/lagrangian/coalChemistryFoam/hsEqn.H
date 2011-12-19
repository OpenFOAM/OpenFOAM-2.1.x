{
    fvScalarMatrix hsEqn
    (
        fvm::ddt(rho, hs)
      + mvConvection->fvmDiv(phi, hs)
      - fvm::laplacian(turbulence->alphaEff(), hs)
     ==
        dpdt
      - (fvc::ddt(rho, K) + fvc::div(phi, K))
      + combustion->Sh()
      + coalParcels.Sh(hs)
      + limestoneParcels.Sh(hs)
      + radiation->Shs(thermo)
      + sources(rho, hs)
    );

    hsEqn.relax();

    sources.constrain(hsEqn);

    hsEqn.solve();

    thermo.correct();

    radiation->correct();

    Info<< "T gas min/max   = " << min(T).value() << ", "
        << max(T).value() << endl;
}
