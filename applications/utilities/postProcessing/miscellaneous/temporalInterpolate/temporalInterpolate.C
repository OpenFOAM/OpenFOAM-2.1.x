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

Description
    Interpolate fields between time-steps e.g. for animation.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"

#include "fvMesh.H"
#include "Time.H"
#include "volMesh.H"
#include "surfaceMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "ReadFields.H"

using namespace Foam;

class fieldInterpolator
{
    Time& runTime_;
    const fvMesh& mesh_;
    const IOobjectList& objects_;
    const HashSet<word>& selectedFields_;
    instant ti_;
    instant ti1_;
    int divisions_;

public:

    fieldInterpolator
    (
        Time& runTime,
        const fvMesh& mesh,
        const IOobjectList& objects,
        const HashSet<word>& selectedFields,
        const instant& ti,
        const instant& ti1,
        int divisions
    )
    :
        runTime_(runTime),
        mesh_(mesh),
        objects_(objects),
        selectedFields_(selectedFields),
        ti_(ti),
        ti1_(ti1),
        divisions_(divisions)
    {}

    template<class GeoFieldType>
    void interpolate();
};


template<class GeoFieldType>
void fieldInterpolator::interpolate()
{
    const word& fieldClassName = GeoFieldType::typeName;

    IOobjectList fields = objects_.lookupClass(fieldClassName);

    if (fields.size())
    {
        Info<< "    " << fieldClassName << "s:";

        forAllConstIter(IOobjectList, fields, fieldIter)
        {
            if
            (
                selectedFields_.empty()
             || selectedFields_.found(fieldIter()->name())
            )
            {
                Info<< " " << fieldIter()->name() << '(';

                GeoFieldType fieldi
                (
                    IOobject
                    (
                        fieldIter()->name(),
                        ti_.name(),
                        fieldIter()->db(),
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    mesh_
                );

                GeoFieldType fieldi1
                (
                    IOobject
                    (
                        fieldIter()->name(),
                        ti1_.name(),
                        fieldIter()->db(),
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    mesh_
                );

                scalar deltaT = (ti1_.value() - ti_.value())/(divisions_ + 1);

                for (int j=0; j<divisions_; j++)
                {
                    instant timej = instant(ti_.value() + (j + 1)*deltaT);

                    runTime_.setTime(timej.name(), 0);

                    Info<< timej.name();

                    if (j < divisions_-1)
                    {
                        Info<< " ";
                    }

                    scalar lambda = scalar(j + 1)/scalar(divisions_ + 1);

                    GeoFieldType fieldj
                    (
                        IOobject
                        (
                            fieldIter()->name(),
                            timej.name(),
                            fieldIter()->db(),
                            IOobject::NO_READ,
                            IOobject::NO_WRITE,
                            false
                        ),
                        (1.0 - lambda)*fieldi + lambda*fieldi1
                    );

                    fieldj.write();
                }

                Info<< ')';
            }
        }

        Info<< endl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:
int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    argList::addOption
    (
        "fields",
        "list",
        "specify a list of fields to be interpolated. Eg, '(U T p)' - "
        "regular expressions not currently supported"
    );
    argList::addOption
    (
        "divisions",
        "integer",
        "specify number of temporal sub-divisions to create (default = 1)."
    );

    #include "setRootCase.H"
    #include "createTime.H"
    runTime.functionObjects().off();

    HashSet<word> selectedFields;
    if (args.optionFound("fields"))
    {
        args.optionLookup("fields")() >> selectedFields;
    }

    int divisions = 1;
    if (args.optionFound("divisions"))
    {
        args.optionLookup("divisions")() >> divisions;
    }

    instantList timeDirs = timeSelector::select0(runTime, args);

    #include "createMesh.H"

    Info<< "Interpolating fields for times:" << endl;

    for (label timei = 0; timei < timeDirs.size() - 1; timei++)
    {
        runTime.setTime(timeDirs[timei], timei);

        // Read objects in time directory
        IOobjectList objects(mesh, runTime.timeName());

        fieldInterpolator interpolator
        (
            runTime,
            mesh,
            objects,
            selectedFields,
            timeDirs[timei],
            timeDirs[timei+1],
            divisions
        );

        // Interpolate vol fields
        interpolator.interpolate<volScalarField>();
        interpolator.interpolate<volVectorField>();
        interpolator.interpolate<volSphericalTensorField>();
        interpolator.interpolate<volSymmTensorField>();
        interpolator.interpolate<volTensorField>();

        // Interpolate surface fields
        interpolator.interpolate<surfaceScalarField>();
        interpolator.interpolate<surfaceVectorField>();
        interpolator.interpolate<surfaceSphericalTensorField>();
        interpolator.interpolate<surfaceSymmTensorField>();
        interpolator.interpolate<surfaceTensorField>();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
