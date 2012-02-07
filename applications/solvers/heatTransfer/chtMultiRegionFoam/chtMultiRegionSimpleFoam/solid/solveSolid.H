{
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix tEqn
        (
            -fvm::laplacian(kappa, T)
        );
        tEqn.relax();
        tEqn.solve();
    }

    Info<< "Min/max T:" << min(T).value() << ' '
        << max(T).value() << endl;
}

thermo.correct();
