/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*----------------------------------------------------------------------------*/

#include "searchableSurfaces.H"
#include "searchableSurfacesQueries.H"
#include "ListOps.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::searchableSurfaces, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct with length.
Foam::searchableSurfaces::searchableSurfaces(const label size)
:
    PtrList<searchableSurface>(size),
    regionNames_(size),
    allSurfaces_(identity(size))
{}


//Foam::searchableSurfaces::searchableSurfaces
//(
//    const IOobject& io,
//    const PtrList<dictionary>& dicts
//)
//:
//    PtrList<searchableSurface>(dicts.size()),
//    regionNames_(dicts.size()),
//    allSurfaces_(identity(dicts.size()))
//{
//    forAll(dicts, surfI)
//    {
//        const dictionary& dict = dicts[surfI];
//
//        // Make IOobject with correct name
//        autoPtr<IOobject> namedIO(io.clone());
//        namedIO().rename(dict.lookup("name"));
//
//        // Create and hook surface
//        set
//        (
//            surfI,
//            searchableSurface::New
//            (
//                dict.lookup("type"),
//                namedIO(),
//                dict
//            )
//        );
//        const searchableSurface& s = operator[](surfI);
//
//        // Construct default region names by prepending surface name
//        // to region name.
//        const wordList& localNames = s.regions();
//
//        wordList globalNames(localNames.size());
//        forAll(localNames, regionI)
//        {
//            globalNames[regionI] = s.name() + '_' + localNames[regionI];
//        }
//
//        // See if dictionary provides any global region names.
//        if (dict.found("regions"))
//        {
//            const dictionary& regionsDict = dict.subDict("regions");
//
//            forAllConstIter(dictionary, regionsDict, iter)
//            {
//                const word& key = iter().keyword();
//
//                if (regionsDict.isDict(key))
//                {
//                    // Get the dictionary for region iter.key()
//                    const dictionary& regionDict = regionsDict.subDict(key);
//
//                    label index = findIndex(localNames, key);
//
//                    if (index == -1)
//                    {
//                        FatalErrorIn
//                        (
//                            "searchableSurfaces::searchableSurfaces"
//                            "( const IOobject&, const dictionary&)"
//                        )   << "Unknown region name " << key
//                            << " for surface " << s.name() << endl
//                            << "Valid region names are " << localNames
//                            << exit(FatalError);
//                    }
//
//                    globalNames[index] = word(regionDict.lookup("name"));
//                }
//            }
//        }
//
//        // Now globalNames contains the names of the regions.
//        Info<< "Surface:" << s.name() << " has regions:"
//            << endl;
//        forAll(globalNames, regionI)
//        {
//            Info<< "    " << globalNames[regionI] << endl;
//        }
//
//        // Create reverse lookup
//        forAll(globalNames, regionI)
//        {
//            regionNames_.insert
//            (
//                globalNames[regionI],
//                labelPair(surfI, regionI)
//            );
//        }
//    }
//}


Foam::searchableSurfaces::searchableSurfaces
(
    const IOobject& io,
    const dictionary& topDict
)
:
    PtrList<searchableSurface>(topDict.size()),
    names_(topDict.size()),
    regionNames_(topDict.size()),
    allSurfaces_(identity(topDict.size()))
{
    label surfI = 0;
    forAllConstIter(dictionary, topDict, iter)
    {
        const word& key = iter().keyword();

        if (!topDict.isDict(key))
        {
            FatalErrorIn
            (
                "searchableSurfaces::searchableSurfaces"
                "( const IOobject&, const dictionary&)"
            )   << "Found non-dictionary entry " << iter()
                << " in top-level dictionary " << topDict
                << exit(FatalError);
        }

        const dictionary& dict = topDict.subDict(key);

        names_[surfI] = key;
        dict.readIfPresent("name", names_[surfI]);

        // Make IOobject with correct name
        autoPtr<IOobject> namedIO(io.clone());
        // Note: we would like to e.g. register triSurface 'sphere.stl' as
        // 'sphere'. Unfortunately
        // no support for having object read from different location than
        // their object name. Maybe have stlTriSurfaceMesh which appends .stl
        // when reading/writing?
        namedIO().rename(key);  // names_[surfI]

        // Create and hook surface
        set
        (
            surfI,
            searchableSurface::New
            (
                dict.lookup("type"),
                namedIO(),
                dict
            )
        );
        const searchableSurface& s = operator[](surfI);

        // Construct default region names by prepending surface name
        // to region name.
        const wordList& localNames = s.regions();

        wordList& rNames = regionNames_[surfI];
        rNames.setSize(localNames.size());

        forAll(localNames, regionI)
        {
            rNames[regionI] = names_[surfI] + '_' + localNames[regionI];
        }

        // See if dictionary provides any global region names.
        if (dict.found("regions"))
        {
            const dictionary& regionsDict = dict.subDict("regions");

            forAllConstIter(dictionary, regionsDict, iter)
            {
                const word& key = iter().keyword();

                if (regionsDict.isDict(key))
                {
                    // Get the dictionary for region iter.keyword()
                    const dictionary& regionDict = regionsDict.subDict(key);

                    label index = findIndex(localNames, key);

                    if (index == -1)
                    {
                        FatalErrorIn
                        (
                            "searchableSurfaces::searchableSurfaces"
                            "( const IOobject&, const dictionary&)"
                        )   << "Unknown region name " << key
                            << " for surface " << s.name() << endl
                            << "Valid region names are " << localNames
                            << exit(FatalError);
                    }

                    rNames[index] = word(regionDict.lookup("name"));
                }
            }
        }

        surfI++;
    }

    // Trim (not really necessary since we don't allow non-dictionary entries)
    PtrList<searchableSurface>::setSize(surfI);
    names_.setSize(surfI);
    regionNames_.setSize(surfI);
    allSurfaces_.setSize(surfI);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::searchableSurfaces::findSurfaceID
(
    const word& wantedName
) const
{
    return findIndex(names_, wantedName);
}


// Find any intersection
void Foam::searchableSurfaces::findAnyIntersection
(
    const pointField& start,
    const pointField& end,
    labelList& hitSurfaces,
    List<pointIndexHit>& hitInfo
) const
{
    searchableSurfacesQueries::findAnyIntersection
    (
        *this,
        allSurfaces_,
        start,
        end,
        hitSurfaces,
        hitInfo
    );
}


//- Find all intersections in order from start to end. Returns for
//  every hit the surface and the hit info.
void Foam::searchableSurfaces::findAllIntersections
(
    const pointField& start,
    const pointField& end,
    labelListList& hitSurfaces,
    List<List<pointIndexHit> >& hitInfo
) const
{
    searchableSurfacesQueries::findAllIntersections
    (
        *this,
        allSurfaces_,
        start,
        end,
        hitSurfaces,
        hitInfo
    );
}


//Find intersections of edge nearest to both endpoints.
void Foam::searchableSurfaces::findNearestIntersection
(
    const pointField& start,
    const pointField& end,
    labelList& surface1,
    List<pointIndexHit>& hit1,
    labelList& surface2,
    List<pointIndexHit>& hit2
) const
{
    searchableSurfacesQueries::findNearestIntersection
    (
        *this,
        allSurfaces_,
        start,
        end,
        surface1,
        hit1,
        surface2,
        hit2
    );
}


// Find nearest. Return -1 or nearest point
void Foam::searchableSurfaces::findNearest
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    labelList& nearestSurfaces,
    List<pointIndexHit>& nearestInfo
) const
{
    searchableSurfacesQueries::findNearest
    (
        *this,
        allSurfaces_,
        samples,
        nearestDistSqr,
        nearestSurfaces,
        nearestInfo
    );
}


//- Calculate bounding box
Foam::boundBox Foam::searchableSurfaces::bounds() const
{
    return searchableSurfacesQueries::bounds
    (
        *this,
        allSurfaces_
    );
}


//- Calculate point which is on a set of surfaces.
Foam::pointIndexHit Foam::searchableSurfaces::facesIntersection
(
    const scalar initDistSqr,
    const scalar convergenceDistSqr,
    const point& start
) const
{
    return searchableSurfacesQueries::facesIntersection
    (
        *this,
        allSurfaces_,
        initDistSqr,
        convergenceDistSqr,
        start
    );
}

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

const Foam::searchableSurface& Foam::searchableSurfaces::operator[]
(
    const word& surfName
) const
{
    const label surfI = findSurfaceID(surfName);

    if (surfI < 0)
    {
        FatalErrorIn
        (
            "searchableSurfaces::operator[](const word&) const"
        )   << "Surface named " << surfName << " not found." << nl
            << "Available surface names: " << names_ << endl
            << abort(FatalError);
    }

    return operator[](surfI);
}


Foam::searchableSurface& Foam::searchableSurfaces::operator[]
(
    const word& surfName
)
{
    const label surfI = findSurfaceID(surfName);

    if (surfI < 0)
    {
        FatalErrorIn
        (
            "searchableSurfaces::operator[](const word&)"
        )   << "Surface named " << surfName << " not found." << nl
            << "Available surface names: " << names_ << endl
            << abort(FatalError);
    }

    return operator[](surfI);
}


// ************************************************************************* //
