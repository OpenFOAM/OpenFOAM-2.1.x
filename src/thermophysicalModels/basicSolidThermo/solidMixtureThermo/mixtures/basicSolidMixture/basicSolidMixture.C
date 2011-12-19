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

#include "basicSolidMixture.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicSolidMixture::basicSolidMixture
(
    const wordList& solidNames,
    const fvMesh& mesh
)
:
    components_(solidNames),
    Y_(components_.size())
{
    forAll(components_, i)
    {
        IOobject header
        (
            "Y" + components_[i],
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ
        );

        // check if field exists and can be read
        if (header.headerOk())
        {
            Y_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "Y" + components_[i],
                        mesh.time().timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh
                )
            );
        }
        else
        {
            volScalarField Ydefault
            (
                IOobject
                (
                    "Ydefault",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );

            Y_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "Y" + components_[i],
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    Ydefault
                )
            );
        }
    }
}


// ************************************************************************* //
