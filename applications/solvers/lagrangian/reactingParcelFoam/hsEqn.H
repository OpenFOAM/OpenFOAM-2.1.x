{
    fvScalarMatrix hEqn
    (
        fvm::ddt(rho, hs)
      + mvConvection->fvmDiv(phi, hs)
      - fvm::laplacian(turbulence->alphaEff(), hs)
     ==
        dpdt
      - (fvc::ddt(rho, K) + fvc::div(phi, K))
      + parcels.Sh(hs)
      + radiation->Shs(thermo)
      + combustion->Sh()
    );

    hEqn.relax();

    hEqn.solve();

    thermo.correct();

    radiation->correct();

    Info<< "T gas min/max   = " << min(T).value() << ", "
        << max(T).value() << endl;
}
