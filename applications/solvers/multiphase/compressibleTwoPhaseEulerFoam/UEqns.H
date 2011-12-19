fvVectorMatrix U1Eqn(U1, U1.dimensions()*dimVol/dimTime);
fvVectorMatrix U2Eqn(U2, U2.dimensions()*dimVol/dimTime);

{
    {
        volTensorField gradU1T(fvc::grad(U1)().T());

        if (kineticTheory.on())
        {
            kineticTheory.solve(gradU1T);
            nuEff1 = kineticTheory.mua()/rho1;
        }
        else // If not using kinetic theory is using Ct model
        {
            nuEff1 = sqr(Ct)*nut2 + nu1;
        }

        volTensorField Rc1
        (
            "Rc1",
            ((2.0/3.0)*I)*(sqr(Ct)*k + nuEff1*tr(gradU1T)) - nuEff1*gradU1T
        );

        if (kineticTheory.on())
        {
            Rc1 -= ((kineticTheory.lambda()/rho1)*tr(gradU1T))*tensor(I);
        }

        U1Eqn =
        (
            (scalar(1) + Cvm*rho2*alpha2/rho1)*
            (
                fvm::ddt(alpha1, U1)
              + fvm::div(alphaPhi1, U1)
              - fvm::Sp(fvc::ddt(alpha1) + fvc::div(alphaPhi1), U1)
            )
          - fvm::laplacian(alpha1*nuEff1, U1)
          + fvc::div(alpha1*Rc1)
         ==
          - fvm::Sp(dragCoeff/rho1, U1)
          - alpha1*alpha2/rho1*(liftForce - Cvm*rho2*DDtU2)
        );

        U1Eqn.relax();
    }

    {
        volTensorField gradU2T(fvc::grad(U2)().T());
        volTensorField Rc2
        (
            "Rc2",
            ((2.0/3.0)*I)*(k + nuEff2*tr(gradU2T)) - nuEff2*gradU2T
        );

        U2Eqn =
        (
            (scalar(1) + Cvm*rho2*alpha1/rho2)*
            (
                fvm::ddt(alpha2, U2)
              + fvm::div(alphaPhi2, U2)
              - fvm::Sp(fvc::ddt(alpha2) + fvc::div(alphaPhi2), U2)
            )
          - fvm::laplacian(alpha2*nuEff2, U2)
          + fvc::div(alpha2*Rc2)
         ==
          - fvm::Sp(dragCoeff/rho2, U2)
          + alpha1*alpha2/rho2*(liftForce + Cvm*rho2*DDtU1)
        );

        U2Eqn.relax();
    }
}
