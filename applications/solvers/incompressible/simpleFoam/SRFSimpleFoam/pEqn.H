{
    p.boundaryField().updateCoeffs();

    volScalarField rAUrel(1.0/UrelEqn().A());
    Urel = rAUrel*UrelEqn().H();
    UrelEqn.clear();

    phi = fvc::interpolate(Urel, "interpolate(HbyA)") & mesh.Sf();
    adjustPhi(phi, Urel, p);

    // Non-orthogonal pressure corrector loop
    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAUrel, p) == fvc::div(phi)
        );

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve();

        if (simple.finalNonOrthogonalIter())
        {
            phi -= pEqn.flux();
        }
    }

    #include "continuityErrs.H"

    // Explicitly relax pressure for momentum corrector
    p.relax();

    // Momentum corrector
    Urel -= rAUrel*fvc::grad(p);
    Urel.correctBoundaryConditions();
    sources.correct(Urel);
}
