/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "isotropicKSolidThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(isotropicKSolidThermo, 0);
    addToRunTimeSelectionTable
    (
        basicSolidThermo,
        isotropicKSolidThermo,
        mesh
    );

    addToRunTimeSelectionTable
    (
        basicSolidThermo,
        isotropicKSolidThermo,
        dictionary
    );

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::isotropicKSolidThermo::isotropicKSolidThermo
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    interpolatedSolidThermo(mesh, typeName + "Coeffs", dict),
    K_
    (
        IOobject
        (
            "K",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimEnergy/dimTime/(dimLength*dimTemperature)
    ),
    KValues_ (Field<scalar>(subDict(typeName + "Coeffs").lookup("KValues")))
{
    correct();
}


Foam::isotropicKSolidThermo::isotropicKSolidThermo(const fvMesh& mesh)
:
    interpolatedSolidThermo(mesh, typeName + "Coeffs"),
    K_
    (
        IOobject
        (
            "K",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimEnergy/dimTime/(dimLength*dimTemperature)
    ),
    KValues_ (Field<scalar>(subDict(typeName + "Coeffs").lookup("KValues")))
{
    correct();
}


void Foam::isotropicKSolidThermo::correct()
{
    // Correct K
    K_.internalField() = interpolateXY
    (
        T_.internalField(),
        TValues_,
        KValues_
    );

    forAll(K_.boundaryField(), patchI)
    {
        K_.boundaryField()[patchI] == interpolateXY
        (
            T_.boundaryField()[patchI],
            TValues_,
            KValues_
        );
    }

    interpolatedSolidThermo::calculate();
}


Foam::tmp<Foam::volScalarField> Foam::isotropicKSolidThermo::K() const
{
    return K_;
}


Foam::tmp<Foam::scalarField> Foam::isotropicKSolidThermo::K
(
    const label patchI
) const
{
    return K_.boundaryField()[patchI];
}


bool Foam::isotropicKSolidThermo::read()
{
    KValues_  = Field<scalar>(subDict(typeName + "Coeffs").lookup("KValues"));
    return true;
}


bool Foam::isotropicKSolidThermo::writeData(Ostream& os) const
{
    os.writeKeyword("KValues") << KValues_ << token::END_STATEMENT << nl;
    bool ok = interpolatedSolidThermo::writeData(os);

    return ok && os.good();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::isotropicKSolidThermo::~isotropicKSolidThermo()
{}


// ************************************************************************* //
