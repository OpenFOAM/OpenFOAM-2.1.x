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

\*---------------------------------------------------------------------------*/

#include "basicSolidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::basicSolidThermo> Foam::basicSolidThermo::New
(
    const fvMesh& mesh
)
{
    if (debug)
    {
        Info<< "basicSolidThermo::New(const fvMesh&): "
            << "constructing basicSolidThermo"
            << endl;
    }

    const word thermoType
    (
        IOdictionary
        (
            IOobject
            (
                "solidThermophysicalProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        ).lookup("thermoType")
    );

    meshConstructorTable::iterator cstrIter =
        meshConstructorTablePtr_->find(thermoType);

    if (cstrIter == meshConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "basicSolidThermo::New(const fvMesh&)"
        )   << "Unknown solidThermo type " << thermoType
            << endl << endl
            << "Valid solidThermo types are :" << endl
            << meshConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<basicSolidThermo>(cstrIter()(mesh));
}


Foam::autoPtr<Foam::basicSolidThermo> Foam::basicSolidThermo::New
(
    const fvMesh& mesh, const dictionary& dict
)
{
    if (debug)
    {
        Info<< "basicSolidThermo::New(const fvMesh&, const dictionary&): "
            << "constructing basicSolidThermo"
            << endl;
    }

    const word thermoType = dict.lookup("thermoType");

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(thermoType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "basicSolidThermo::New(const fvMesh&, const dictionary&)"
        )   << "Unknown solidThermo type " << thermoType
            << endl << endl
            << "Valid solidThermo types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<basicSolidThermo>(cstrIter()(mesh, dict));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
