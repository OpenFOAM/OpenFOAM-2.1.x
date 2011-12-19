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
    patchIntegrate

Description
    Calculates the integral of the specified field over the specified patch.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
#   include "addRegionOption.H"
    timeSelector::addOptions();
    argList::validArgs.append("fieldName");
    argList::validArgs.append("patchName");
#   include "setRootCase.H"
#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createNamedMesh.H"

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

            // Give patch area
            Info<< "    Area vector of patch "
                << patchName << '[' << patchI << ']' << " = "
                << gSum(mesh.Sf().boundaryField()[patchI]) << endl;
            Info<< "    Area magnitude of patch "
                << patchName << '[' << patchI << ']' << " = "
                << gSum(mesh.magSf().boundaryField()[patchI]) << endl;

            // Read field and calc integral
            if (fieldHeader.headerClassName() == volScalarField::typeName)
            {
                Info<< "    Reading " << volScalarField::typeName << " "
                    << fieldName << endl;

                volScalarField field(fieldHeader, mesh);

                Info<< "    Integral of " << fieldName
                    << " over vector area of patch "
                    << patchName << '[' << patchI << ']' << " = "
                    << gSum
                       (
                           mesh.Sf().boundaryField()[patchI]
                          *field.boundaryField()[patchI]
                       )
                    << nl;

                Info<< "    Integral of " << fieldName
                    << " over area magnitude of patch "
                    << patchName << '[' << patchI << ']' << " = "
                    << gSum
                       (
                           mesh.magSf().boundaryField()[patchI]
                          *field.boundaryField()[patchI]
                       )
                    << nl;
            }
            else if
            (
                fieldHeader.headerClassName() == surfaceScalarField::typeName
            )
            {
                Info<< "    Reading " << surfaceScalarField::typeName << " "
                    << fieldName << endl;

                surfaceScalarField field(fieldHeader, mesh);
                scalar sumField = gSum(field.boundaryField()[patchI]);

                Info<< "    Integral of " << fieldName << " over patch "
                    << patchName << '[' << patchI << ']' << " = "
                    << sumField << nl;
            }
            else
            {
                FatalError
                    << "Only possible to integrate "
                    << volScalarField::typeName << "s "
                    << "and " << surfaceScalarField::typeName << "s"
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
