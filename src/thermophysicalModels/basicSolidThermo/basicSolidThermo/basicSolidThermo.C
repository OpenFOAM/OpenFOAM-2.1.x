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

#include "basicSolidThermo.H"
#include "fvMesh.H"
#include "HashTable.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(basicSolidThermo, 0);
    defineRunTimeSelectionTable(basicSolidThermo, mesh);
    defineRunTimeSelectionTable(basicSolidThermo, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicSolidThermo::basicSolidThermo(const fvMesh& mesh)
:
    IOdictionary
    (
        IOobject
        (
            "solidThermophysicalProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    T_
    (
        IOobject
        (
            "T",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    rho_
    (
        IOobject
        (
            "rho",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimMass/dimVolume
    ),
    kappa_
    (
        IOobject
        (
            "kappa",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimless/dimLength
    ),
    sigmaS_
    (
        IOobject
        (
            "sigmaS",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimless/dimLength
    ),
    emissivity_
    (
        IOobject
        (
            "emissivity",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimless
    )

{}


Foam::basicSolidThermo::basicSolidThermo
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    IOdictionary
    (
        IOobject
        (
            "solidThermophysicalProperties",
            mesh.time().constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dict
    ),
    mesh_(mesh),
    T_
    (
        IOobject
        (
            "T",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    rho_
    (
        IOobject
        (
            "rho",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimMass/dimVolume
    ),
    kappa_
    (
        IOobject
        (
            "kappa",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimless/dimLength
    ),
    sigmaS_
    (
        IOobject
        (
            "sigmaS",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimless/dimLength
    ),
    emissivity_
    (
        IOobject
        (
            "emissivity",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimless
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::basicSolidThermo::~basicSolidThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::volScalarField& Foam::basicSolidThermo::T()
{
    return T_;
}


const Foam::volScalarField& Foam::basicSolidThermo::T() const
{
    return T_;
}


const Foam::volScalarField& Foam::basicSolidThermo::rho() const
{
    return rho_;
}


Foam::volScalarField& Foam::basicSolidThermo::rho()
{
    return rho_;
}


const Foam::volScalarField& Foam::basicSolidThermo::kappa() const
{
    return kappa_;
}


const Foam::volScalarField& Foam::basicSolidThermo::sigmaS() const
{
    return sigmaS_;
}


const Foam::volScalarField& Foam::basicSolidThermo::emissivity() const
{
    return emissivity_;
}


Foam::basicSolidMixture& Foam::basicSolidThermo::composition()
{
    notImplemented("basicSolidThermo::composition()");
    return *reinterpret_cast<basicSolidMixture*>(0);
}


const Foam::basicSolidMixture& Foam::basicSolidThermo::composition() const
{
    notImplemented("basicSolidThermo::composition() const");
    return *reinterpret_cast<const basicSolidMixture*>(0);
}


Foam::tmp<Foam::volScalarField> Foam::basicSolidThermo::hs() const
{
    notImplemented("basicSolidThermo::hs() const");
    return volScalarField::null();
}


Foam::tmp<Foam::scalarField> Foam::basicSolidThermo::hs(const label patchI)
const
{
    notImplemented("basicSolidThermo::hs(const label) const");
    return scalarField::null();
}


bool Foam::basicSolidThermo::read()
{
    return regIOobject::read();
}


bool Foam::basicSolidThermo::writeData(Ostream& os) const
{
    return true;
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const basicSolidThermo& s)
{
    s.writeData(os);
    return os;
}


// ************************************************************************* //
