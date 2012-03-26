{
    // rho1 = rho10 + psi1*p;
    // rho2 = rho20 + psi2*p;

    // tmp<fvScalarMatrix> pEqnComp1;
    // tmp<fvScalarMatrix> pEqnComp2;

    // //if (transonic)
    // //{
    // //}
    // //else
    // {
    //     surfaceScalarField phid1("phid1", fvc::interpolate(psi1)*phi1);
    //     surfaceScalarField phid2("phid2", fvc::interpolate(psi2)*phi2);

    //     pEqnComp1 =
    //         fvc::ddt(rho1) + psi1*correction(fvm::ddt(p))
    //       + fvc::div(phid1, p)
    //       - fvc::Sp(fvc::div(phid1), p);

    //     pEqnComp2 =
    //         fvc::ddt(rho2) + psi2*correction(fvm::ddt(p))
    //       + fvc::div(phid2, p)
    //       - fvc::Sp(fvc::div(phid2), p);
    // }

    PtrList<surfaceScalarField> alphafs(fluid.phases().size());
    PtrList<volVectorField> U0s(fluid.phases().size());
    PtrList<surfaceScalarField> phi0s(fluid.phases().size());
    PtrList<volScalarField> rAUs(fluid.phases().size());
    PtrList<surfaceScalarField> rAlphaAUfs(fluid.phases().size());

    phasei = 0;
    forAllIter(PtrDictionary<phaseModel>, fluid.phases(), iter)
    {
        phaseModel& phase = iter();

        U0s.set(phasei, new volVectorField(phase.U()));
        phi0s.set(phasei, new surfaceScalarField(phase.phi()));
        mrfZones.absoluteFlux(phi0s[phasei]);

        phasei++;

    }

    surfaceScalarField phi0
    (
        IOobject
        (
            "phiHbyA",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("phiHbyA", dimArea*dimVelocity, 0)
    );

    phasei = 0;
    forAllIter(PtrDictionary<phaseModel>, fluid.phases(), iter)
    {
        phaseModel& phase = iter();
        const volScalarField& alpha = phase;

        alphafs.set(phasei, fvc::interpolate(alpha).ptr());
        rAUs.set(phasei, (1.0/UEqns[phasei].A()).ptr());
        rAlphaAUfs.set(phasei, fvc::interpolate(alpha*rAUs[phasei]).ptr());

        phase.U() = rAUs[phasei]*UEqns[phasei].H();

        phase.phi().oldTime();
        mrfZones.absoluteFlux(phase.phi().oldTime());
        mrfZones.absoluteFlux(phase.phi());
        phase.phi() =
        (
            (fvc::interpolate(phase.U()) & mesh.Sf())
          + fvc::ddtPhiCorr(rAUs[phasei], alpha, phase.U(), phase.phi())
        );
        mrfZones.relativeFlux(phase.phi());
        mrfZones.relativeFlux(phase.phi().oldTime());
        surfaceScalarField pphi0("pphi0", 1.0*phase.phi());
        pphi0 += rAlphaAUfs[phasei]*(g & mesh.Sf());

        multiphaseSystem::dragModelTable::const_iterator dmIter =
            fluid.dragModels().begin();
        multiphaseSystem::dragCoeffFields::const_iterator dcIter =
            dragCoeffs().begin();
        for
        (
            ;
            dmIter != fluid.dragModels().end() && dcIter != dragCoeffs().end();
            ++dmIter, ++dcIter
        )
        {
            if
            (
                &phase == &dmIter()->phase1()
             || &phase == &dmIter()->phase2()
            )
            {
                int phasej = 0;
                forAllIter(PtrDictionary<phaseModel>, fluid.phases(), iter2)
                {
                    phaseModel& phase2 = iter2();

                    if
                    (
                        (
                            &phase == &dmIter()->phase1()
                         && &phase2 == &dmIter()->phase2()
                        )
                     || (
                            &phase == &dmIter()->phase2()
                         && &phase2 == &dmIter()->phase1()
                        )
                    ) break;

                    phasej++;
                }

                pphi0 +=
                    fvc::interpolate
                    ((1.0/phase.rho())*rAUs[phasei]*(*dcIter()))
                   *phi0s[phasej];
            }
        }

        phi0 += alphafs[phasei]*pphi0;
        phase.phi() = pphi0;

        phasei++;
    }

    phi = dimensionedScalar("phi", phi.dimensions(), 0);
    surfaceScalarField Dp
    (
        IOobject
        (
            "Dp",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("Dp", dimensionSet(-1, 3, 1, 0, 0), 0)
    );

    phasei = 0;
    forAllIter(PtrDictionary<phaseModel>, fluid.phases(), iter)
    {
        phaseModel& phase = iter();
        phi += alphafs[phasei]*phase.phi();
        Dp += alphafs[phasei]*rAlphaAUfs[phasei]/phase.rho();

        phasei++;
    }
    Dp = mag(Dp);
    adjustPhi(phi0, U, p);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqnIncomp
        (
            fvc::div(phi0)
          - fvm::laplacian(Dp, p)
        );

        solve
        (
            // (
            //     (alpha1/rho1)*pEqnComp1()
            //   + (alpha2/rho2)*pEqnComp2()
            // ) +
            pEqnIncomp,
            mesh.solver(p.select(pimple.finalInnerIter()))
        );

        if (pimple.finalNonOrthogonalIter())
        {
            surfaceScalarField mSfGradp(pEqnIncomp.flux()/Dp);

            phasei = 0;
            phi = dimensionedScalar("phi", phi.dimensions(), 0);
            forAllIter(PtrDictionary<phaseModel>, fluid.phases(), iter)
            {
                phaseModel& phase = iter();

                phase.phi() += rAlphaAUfs[phasei]*mSfGradp/phase.rho();
                phi += alphafs[phasei]*phase.phi();

                phasei++;
            }

            // dgdt =

            // (
            //     pos(alpha2)*(pEqnComp2 & p)/rho2
            //   - pos(alpha1)*(pEqnComp1 & p)/rho1
            // );

            p.relax();
            mSfGradp = pEqnIncomp.flux()/Dp;

            U = dimensionedVector("U", dimVelocity, vector::zero);

            phasei = 0;
            forAllIter(PtrDictionary<phaseModel>, fluid.phases(), iter)
            {
                phaseModel& phase = iter();
                const volScalarField& alpha = phase;

                phase.U() +=
                    fvc::reconstruct
                    (
                        rAlphaAUfs[phasei]*(g & mesh.Sf())
                      + rAlphaAUfs[phasei]*mSfGradp/phase.rho()
                    );

                multiphaseSystem::dragModelTable::const_iterator dmIter =
                    fluid.dragModels().begin();
                multiphaseSystem::dragCoeffFields::const_iterator dcIter =
                    dragCoeffs().begin();
                for
                (
                    ;
                    dmIter != fluid.dragModels().end()
                 && dcIter != dragCoeffs().end();
                    ++dmIter, ++dcIter
                )
                {
                    if
                    (
                        &phase == &dmIter()->phase1()
                     || &phase == &dmIter()->phase2()
                    )
                    {
                        int phasej = 0;
                        forAllIter
                        (
                            PtrDictionary<phaseModel>,
                            fluid.phases(),
                            iter2
                        )
                        {
                            phaseModel& phase2 = iter2();

                            if
                            (
                                (
                                    &phase == &dmIter()->phase1()
                                 && &phase2 == &dmIter()->phase2()
                                )
                             || (
                                    &phase == &dmIter()->phase2()
                                 && &phase2 == &dmIter()->phase1()
                                 )
                            ) break;

                            phasej++;
                        }

                        // phase.U() +=
                        //     (1.0/phase.rho())*rAUs[phasei]*(*dcIter())
                        //    *U0s[phasej];

                        phase.U() += fvc::reconstruct
                        (
                            fvc::interpolate
                            (
                                (1.0/phase.rho())*rAUs[phasei]*(*dcIter())
                            )*phi0s[phasej]
                        );
                    }
                }

                phase.U().correctBoundaryConditions();

                U += alpha*phase.U();

                phasei++;
            }
        }
    }

    //p = max(p, pMin);

    #include "continuityErrs.H"

    // rho1 = rho10 + psi1*p;
    // rho2 = rho20 + psi2*p;

    // Dp1Dt = fvc::DDt(phi1, p);
    // Dp2Dt = fvc::DDt(phi2, p);
}
