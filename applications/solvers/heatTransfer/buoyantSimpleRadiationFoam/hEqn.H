{
    fvScalarMatrix hEqn
    (
        fvm::div(phi, h)
      - fvm::Sp(fvc::div(phi), h)
      - fvm::laplacian(turbulence->alphaEff(), h)
     ==
      - fvc::div(phi, 0.5*magSqr(U), "div(phi,K)")
      + radiation->Sh(thermo)
    );

    hEqn.relax();

    hEqn.solve();

    thermo.correct();

    radiation->correct();
}
