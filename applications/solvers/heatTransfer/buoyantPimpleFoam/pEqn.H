{
    rho = thermo.rho();

    // Thermodynamic density needs to be updated by psi*d(p) after the
    // pressure solution - done in 2 parts. Part 1:
    thermo.rho() -= psi*p_rgh;

    volScalarField rAU(1.0/UEqn.A());
    surfaceScalarField rhorAUf("(rho*(1|A(U)))", fvc::interpolate(rho*rAU));

    U = rAU*UEqn.H();

    phi = fvc::interpolate(rho)*
    (
        (fvc::interpolate(U) & mesh.Sf())
      + fvc::ddtPhiCorr(rAU, rho, U, phi)
    );

    surfaceScalarField buoyancyPhi(-rhorAUf*ghf*fvc::snGrad(rho)*mesh.magSf());
    phi += buoyancyPhi;

    fvScalarMatrix p_rghDDtEqn
    (
        fvc::ddt(rho) + psi*correction(fvm::ddt(p_rgh))
      + fvc::div(phi)
    );

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            p_rghDDtEqn
          - fvm::laplacian(rhorAUf, p_rgh)
        );

        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            // Calculate the conservative fluxes
            phi += p_rghEqn.flux();

            // Explicitly relax pressure for momentum corrector
            p_rgh.relax();

            // Correct the momentum source with the pressure gradient flux
            // calculated from the relaxed pressure
            U += rAU*fvc::reconstruct((buoyancyPhi + p_rghEqn.flux())/rhorAUf);
            U.correctBoundaryConditions();
            K = 0.5*magSqr(U);
        }
    }

    p = p_rgh + rho*gh;

    // Second part of thermodynamic density update
    thermo.rho() += psi*p_rgh;

    dpdt = fvc::ddt(p);

    #include "rhoEqn.H"
    #include "compressibleContinuityErrs.H"
}
