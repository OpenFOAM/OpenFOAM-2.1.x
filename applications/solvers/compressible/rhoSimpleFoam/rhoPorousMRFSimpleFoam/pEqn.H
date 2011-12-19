if (pressureImplicitPorosity)
{
    U = trTU()&UEqn().H();
}
else
{
    U = trAU()*UEqn().H();
}

UEqn.clear();

bool closedVolume = false;

if (simple.transonic())
{
    surfaceScalarField phid
    (
        "phid",
        fvc::interpolate(psi)*(fvc::interpolate(U) & mesh.Sf())
    );
    mrfZones.relativeFlux(fvc::interpolate(psi), phid);

    while (simple.correctNonOrthogonal())
    {
        tmp<fvScalarMatrix> tpEqn;

        if (pressureImplicitPorosity)
        {
            tpEqn = (fvc::div(phid, p) - fvm::laplacian(rho*trTU(), p));
        }
        else
        {
            tpEqn = (fvc::div(phid, p) - fvm::laplacian(rho*trAU(), p));
        }

        tpEqn().setReference(pRefCell, pRefValue);

        tpEqn().solve();

        if (simple.finalNonOrthogonalIter())
        {
            phi == tpEqn().flux();
        }
    }
}
else
{
    phi = fvc::interpolate(rho*U) & mesh.Sf();
    mrfZones.relativeFlux(fvc::interpolate(rho), phi);

    closedVolume = adjustPhi(phi, U, p);

    while (simple.correctNonOrthogonal())
    {
        tmp<fvScalarMatrix> tpEqn;

        if (pressureImplicitPorosity)
        {
            tpEqn = (fvm::laplacian(rho*trTU(), p) == fvc::div(phi));
        }
        else
        {
            tpEqn = (fvm::laplacian(rho*trAU(), p) == fvc::div(phi));
        }

        tpEqn().setReference(pRefCell, pRefValue);

        tpEqn().solve();

        if (simple.finalNonOrthogonalIter())
        {
            phi -= tpEqn().flux();
        }
    }
}

#include "incompressible/continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

if (pressureImplicitPorosity)
{
    U -= trTU()&fvc::grad(p);
}
else
{
    U -= trAU()*fvc::grad(p);
}

U.correctBoundaryConditions();

// For closed-volume cases adjust the pressure and density levels
// to obey overall mass continuity
if (closedVolume)
{
    p += (initialMass - fvc::domainIntegrate(psi*p))
        /fvc::domainIntegrate(psi);
}

rho = thermo.rho();
rho = max(rho, rhoMin);
rho = min(rho, rhoMax);
rho.relax();
Info<< "rho max/min : " << max(rho).value() << " " << min(rho).value() << endl;
