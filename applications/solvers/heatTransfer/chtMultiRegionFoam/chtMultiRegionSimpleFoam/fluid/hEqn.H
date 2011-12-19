{
    fvScalarMatrix hEqn
    (
        fvm::div(phi, h)
      - fvm::Sp(fvc::div(phi), h)
      - fvm::laplacian(turb.alphaEff(), h)
     ==
      - fvc::div(phi, 0.5*magSqr(U), "div(phi,K)")
      + rad.Sh(thermo)
    );

    hEqn.relax();

    hEqn.solve();

    thermo.correct();

    rad.correct();

    Info<< "Min/max T:" << min(thermo.T()).value() << ' '
        << max(thermo.T()).value() << endl;
}
