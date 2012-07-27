fvVectorMatrix U1Eqn(U1, U1.dimensions()*dimVol/dimTime);
fvVectorMatrix U2Eqn(U2, U2.dimensions()*dimVol/dimTime);

{
    volTensorField Rc1(-nuEff1*(T(fvc::grad(U1))));
    Rc1 = Rc1 + (2.0/3.0)*sqr(Ct)*I*k - (2.0/3.0)*I*tr(Rc1);

    surfaceScalarField phiR1
    (
      - fvc::interpolate(nuEff1)
       *mesh.magSf()*fvc::snGrad(alpha1)
       /fvc::interpolate(alpha1 + scalar(0.001))
    );

    U1Eqn =
    (
        (scalar(1) + Cvm*rho2*alpha2/rho1)*
        (
            fvm::ddt(U1)
          + fvm::div(phi1, U1, "div(phi1,U1)")
          - fvm::Sp(fvc::div(phi1), U1)
        )

      - fvm::laplacian(nuEff1, U1)
      + fvc::div(Rc1)

      + fvm::div(phiR1, U1, "div(phi1,U1)")
      - fvm::Sp(fvc::div(phiR1), U1)
      + (fvc::grad(alpha1)/(fvc::average(alpha1) + scalar(0.001)) & Rc1)
     ==
    //  g                        // Buoyancy term transfered to p-equation
      - fvm::Sp(alpha2/rho1*dragCoef, U1)
    //+ alpha2/rho1*dragCoef*U2    // Explicit drag transfered to p-equation
      - alpha2/rho1*(liftCoeff - Cvm*rho2*DDtU2)
    );

    U1Eqn.relax();


    volTensorField Rc2(-nuEff2*T(fvc::grad(U2)));
    Rc2 = Rc2 + (2.0/3.0)*I*k - (2.0/3.0)*I*tr(Rc2);

    surfaceScalarField phiR2
    (
       - fvc::interpolate(nuEff2)
        *mesh.magSf()*fvc::snGrad(alpha2)
        /fvc::interpolate(alpha2 + scalar(0.001))
    );

    U2Eqn =
    (
        (scalar(1) + Cvm*rho2*alpha1/rho2)*
        (
            fvm::ddt(U2)
          + fvm::div(phi2, U2, "div(phi2,U2)")
          - fvm::Sp(fvc::div(phi2), U2)
        )

      - fvm::laplacian(nuEff2, U2)
      + fvc::div(Rc2)

      + fvm::div(phiR2, U2, "div(phi2,U2)")
      - fvm::Sp(fvc::div(phiR2), U2)

      + (fvc::grad(alpha2)/(fvc::average(alpha2) + scalar(0.001)) & Rc2)
     ==
    //  g                        // Buoyancy term transfered to p-equation
      - fvm::Sp(alpha1/rho2*dragCoef, U2)
    //+ alpha1/rho2*dragCoef*U1   // Explicit drag transfered to p-equation
      + alpha1/rho2*(liftCoeff + Cvm*rho2*DDtU1)
    );

    U2Eqn.relax();
}
