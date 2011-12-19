    MRFZones mrfZones(mesh);
    mrfZones.correctBoundaryVelocity(U);

    porousZones pZones(mesh);
    Switch pressureImplicitPorosity(false);
