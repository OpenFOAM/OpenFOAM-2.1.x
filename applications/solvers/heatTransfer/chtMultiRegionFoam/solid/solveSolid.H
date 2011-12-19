if (finalIter)
{
    mesh.data::add("finalIteration", true);
}

{
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        tmp<fvScalarMatrix> TEqn
        (
            fvm::ddt(rho*cp, T)
          - fvm::laplacian(kappa, T)
        );
        TEqn().relax();
        TEqn().solve(mesh.solver(T.select(finalIter)));
    }

    Info<< "Min/max T:" << min(T) << ' ' << max(T) << endl;
}

thermo.correct();

if (finalIter)
{
    mesh.data::remove("finalIteration");
}
