word dictName(args.optionLookupOrDefault<word>("dict", "particleTrackDict"));

IOdictionary propsDict
(
    IOobject
    (
        dictName,
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED
    )
);

word cloudName(propsDict.lookup("cloudName"));

List<word> userFields(propsDict.lookup("fields"));