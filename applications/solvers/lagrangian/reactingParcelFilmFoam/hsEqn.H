{
    fvScalarMatrix hsEqn
    (
        fvm::ddt(rho, hs)
      + mvConvection->fvmDiv(phi, hs)
      - fvm::laplacian(turbulence->alphaEff(), hs)
     ==
        dpdt
      - (fvc::ddt(rho, K) + fvc::div(phi, K))
      + parcels.Sh(hs)
      + surfaceFilm.Sh()
      + radiation->Shs(thermo)
      + combustion->Sh()
    );

    hsEqn.relax();

    hsEqn.solve();

    thermo.correct();

    radiation->correct();

    Info<< "min/max(T) = " << min(T).value() << ", " << max(T).value() << endl;
}
