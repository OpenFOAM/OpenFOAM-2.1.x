{
    rho = thermo.rho();

    // Thermodynamic density needs to be updated by psi*d(p) after the
    // pressure solution - done in 2 parts. Part 1:
    thermo.rho() -= psi*p;

    volScalarField rAU(1.0/UEqn.A());
    U = rAU*(UEqn == sources(rho, U))().H();

    if (pZones.size() > 0)
    {
        // ddtPhiCorr not well defined for cases with porosity
        phi = fvc::interpolate(rho)*(fvc::interpolate(U) & mesh.Sf());
    }
    else
    {
        phi =
            fvc::interpolate(rho)
           *(
                (fvc::interpolate(U) & mesh.Sf())
              + fvc::ddtPhiCorr(rAU, rho, U, phi)
            );
    }

    fvScalarMatrix pDDtEqn
    (
        fvc::ddt(rho) + psi*correction(fvm::ddt(p))
      + fvc::div(phi)
     ==
        parcels.Srho()
      + sources(psi, p, rho.name())
    );


    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            pDDtEqn
          - fvm::laplacian(rho*rAU, p)
        );

        sources.constrain(pDDtEqn, rho.name());

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi += pEqn.flux();
        }
    }

    // Second part of thermodynamic density update
    thermo.rho() += psi*p;

    #include "rhoEqn.H"
    #include "compressibleContinuityErrs.H"

    U -= rAU*fvc::grad(p);
    U.correctBoundaryConditions();
    sources.correct(U);

    K = 0.5*magSqr(U);

    dpdt = fvc::ddt(p);
}
