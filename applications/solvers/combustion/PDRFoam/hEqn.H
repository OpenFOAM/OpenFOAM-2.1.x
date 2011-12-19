{
    solve
    (
        betav*fvm::ddt(rho, h)
      + mvConvection->fvmDiv(phi, h)
      - fvm::laplacian(Db, h)
     ==
        betav*dpdt
      - betav*(fvc::ddt(rho, K) + fvc::div(phi, K))
    );

    thermo.correct();
}
