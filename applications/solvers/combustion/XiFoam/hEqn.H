{
    fvScalarMatrix hEqn
    (
        fvm::ddt(rho, h)
      + mvConvection->fvmDiv(phi, h)
      - fvm::laplacian(turbulence->alphaEff(), h)
     ==
        dpdt
      - (fvc::ddt(rho, K) + fvc::div(phi, K))
    );

    hEqn.relax();
    hEqn.solve();

    thermo.correct();
}
