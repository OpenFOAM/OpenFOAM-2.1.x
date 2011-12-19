    // Initialise solid field pointer lists
    PtrList<basicSolidThermo> thermos(solidRegions.size());

    // Populate solid field pointer lists
    forAll(solidRegions, i)
    {
        Info<< "*** Reading solid mesh thermophysical properties for region "
            << solidRegions[i].name() << nl << endl;

        Info<< "    Adding to thermos\n" << endl;
        thermos.set(i, basicSolidThermo::New(solidRegions[i]));
    }
