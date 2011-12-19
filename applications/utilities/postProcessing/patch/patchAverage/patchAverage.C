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

Application
    patchAverage

Description
    Calculates the average of the specified field over the specified patch.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    argList::validArgs.append("fieldName");
    argList::validArgs.append("patchName");
#   include "setRootCase.H"
#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createMesh.H"

    const word fieldName = args[1];
    const word patchName = args[2];

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        IOobject fieldHeader
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        // Check field exists
        if (fieldHeader.headerOk())
        {
            mesh.readUpdate();

            const label patchI = mesh.boundaryMesh().findPatchID(patchName);
            if (patchI < 0)
            {
                FatalError
                    << "Unable to find patch " << patchName << nl
                    << exit(FatalError);
            }

            if (fieldHeader.headerClassName() == "volScalarField")
            {
                Info<< "    Reading volScalarField " << fieldName << endl;
                volScalarField field(fieldHeader, mesh);

                scalar area = gSum(mesh.magSf().boundaryField()[patchI]);
                scalar sumField = 0;

                if (area > 0)
                {
                    sumField = gSum
                    (
                        mesh.magSf().boundaryField()[patchI]
                      * field.boundaryField()[patchI]
                    ) / area;
                }

                Info<< "    Average of " << fieldName << " over patch "
                    << patchName << '[' << patchI << ']' << " = "
                    << sumField << endl;
            }
            else
            {
                FatalError
                    << "Only possible to average volScalarFields "
                    << nl << exit(FatalError);
            }
        }
        else
        {
            Info<< "    No field " << fieldName << endl;
        }

        Info<< endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
