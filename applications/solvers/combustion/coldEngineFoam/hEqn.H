{
    solve
    (
        fvm::ddt(rho, h)
      + fvm::div(phi, h)
      - fvm::laplacian(turbulence->alphaEff(), h)
     ==
      - fvc::div(phi, 0.5*magSqr(U))
    );

    thermo.correct();
}
