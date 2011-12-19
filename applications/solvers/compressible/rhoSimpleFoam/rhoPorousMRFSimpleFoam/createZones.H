    MRFZones mrfZones(mesh);
    mrfZones.correctBoundaryVelocity(U);

    thermalPorousZones pZones(mesh);
    Switch pressureImplicitPorosity(false);

    // nUCorrectors used for pressureImplicitPorosity
    int nUCorr = 0;
    if (pZones.size())
    {
        // nUCorrectors for pressureImplicitPorosity
        simple.dict().readIfPresent("nUCorrectors", nUCorr);

        if (nUCorr > 0)
        {
            pressureImplicitPorosity = true;
            Info<< "Using pressure implicit porosity" << endl;
        }
        else
        {
            Info<< "Using pressure explicit porosity" << endl;
        }
    }
