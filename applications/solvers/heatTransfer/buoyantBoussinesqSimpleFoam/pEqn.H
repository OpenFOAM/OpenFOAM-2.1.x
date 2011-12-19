{
    volScalarField rAU("rAU", 1.0/UEqn().A());
    surfaceScalarField rAUf("(1|A(U))", fvc::interpolate(rAU));

    U = rAU*UEqn().H();
    UEqn.clear();

    phi = fvc::interpolate(U) & mesh.Sf();
    adjustPhi(phi, U, p_rgh);

    surfaceScalarField buoyancyPhi(rAUf*ghf*fvc::snGrad(rhok)*mesh.magSf());
    phi -= buoyancyPhi;

    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUf, p_rgh) == fvc::div(phi)
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        p_rghEqn.solve();

        if (simple.finalNonOrthogonalIter())
        {
            // Calculate the conservative fluxes
            phi -= p_rghEqn.flux();

            // Explicitly relax pressure for momentum corrector
            p_rgh.relax();

            // Correct the momentum source with the pressure gradient flux
            // calculated from the relaxed pressure
            U -= rAU*fvc::reconstruct((buoyancyPhi + p_rghEqn.flux())/rAUf);
            U.correctBoundaryConditions();
        }
    }

    #include "continuityErrs.H"

    p = p_rgh + rhok*gh;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rhok*gh;
    }
}
