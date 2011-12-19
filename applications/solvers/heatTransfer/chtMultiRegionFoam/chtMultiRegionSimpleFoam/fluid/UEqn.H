    // Solve the Momentum equation
    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + turb.divDevRhoReff(U)
    );

    UEqn().relax();

    solve
    (
        UEqn()
     ==
        fvc::reconstruct
        (
            (
              - ghf*fvc::snGrad(rho)
              - fvc::snGrad(p_rgh)
            )*mesh.magSf()
        )
    );
