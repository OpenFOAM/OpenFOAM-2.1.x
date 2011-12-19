{
    fvScalarMatrix hsEqn
    (
        fvm::ddt(rho, hs)
      + mvConvection->fvmDiv(phi, hs)
      - fvm::laplacian(turbulence->alphaEff(), hs)
     ==
      - (fvc::ddt(rho, K) + fvc::div(phi, K))
      + parcels.Sh(hs)
      + radiation->Shs(thermo)
      + combustion->Sh()
      + sources(rho, hs)
    );

    if (pressureWorkTimeDerivative)
    {
        hsEqn -= dpdt;
    }

    hsEqn.relax();

    sources.constrain(hsEqn);

    hsEqn.solve();

    thermo.correct();

    radiation->correct();

    Info<< "T gas min/max   = " << min(T).value() << ", "
        << max(T).value() << endl;
}
