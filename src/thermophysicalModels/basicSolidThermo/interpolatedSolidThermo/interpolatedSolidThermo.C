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

#include "interpolatedSolidThermo.H"
#include "interpolateXY.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::interpolatedSolidThermo::rho
(
    const label patchI
) const
{
    return tmp<scalarField>
    (
        new scalarField
        (
            interpolateXY
            (
                T_.boundaryField()[patchI],
                TValues_,
                rhoValues_
            )
        )
    );
}


Foam::tmp<Foam::scalarField> Foam::interpolatedSolidThermo::emissivity
(
    const label patchI
) const
{
    return tmp<scalarField>
    (
        new scalarField
        (
            interpolateXY
            (
                T_.boundaryField()[patchI],
                TValues_,
                emissivityValues_
            )
        )
    );
}


Foam::tmp<Foam::scalarField> Foam::interpolatedSolidThermo::kappa
(
    const label patchI
) const
{
    return tmp<scalarField>
    (
        new scalarField
        (
            interpolateXY
            (
                T_.boundaryField()[patchI],
                TValues_,
                kappaValues_
            )
        )
    );
}


Foam::tmp<Foam::scalarField> Foam::interpolatedSolidThermo::sigmaS
(
    const label patchI
) const
{
    return tmp<scalarField>
    (
        new scalarField
        (
            interpolateXY
            (
                T_.boundaryField()[patchI],
                TValues_,
                sigmaSValues_
            )
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interpolatedSolidThermo::interpolatedSolidThermo
(
    const fvMesh& mesh,
    const word dictName
 )
:
    basicSolidThermo(mesh),
    interpolateSolid(subDict(dictName)),
    dict_(subDict(dictName))
{
    calculate();
}


Foam::interpolatedSolidThermo::interpolatedSolidThermo
(
    const fvMesh& mesh,
    const word dictName,
    const dictionary& dict
 )
:
    basicSolidThermo(mesh, dict),
    interpolateSolid(subDict(dictName)),
    dict_(subDict(dictName))
{
    calculate();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interpolatedSolidThermo::~interpolatedSolidThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::interpolatedSolidThermo::calculate()
{
    // Correct rho
    rho_.internalField() = interpolateXY
    (
        T_.internalField(),
        TValues_,
        rhoValues_
    );

    forAll(rho_.boundaryField(), patchI)
    {
        rho_.boundaryField()[patchI] == this->rho(patchI)();
    }

    // Correct emissivity
    emissivity_.internalField() = interpolateXY
    (
        T_.internalField(),
        TValues_,
        emissivityValues_
    );

    forAll(emissivity_.boundaryField(), patchI)
    {
        emissivity_.boundaryField()[patchI] == this->emissivity(patchI)();
    }


    // Correct absorptivity
    kappa_.internalField() = interpolateXY
    (
        T_.internalField(),
        TValues_,
        kappaValues_
    );

    forAll(kappa_.boundaryField(), patchI)
    {
        kappa_.boundaryField()[patchI] == this->kappa(patchI)();
    }


    // Correct scatter
    sigmaS_.internalField() = interpolateXY
    (
        T_.internalField(),
        TValues_,
        sigmaSValues_
    );

    forAll(sigmaS_.boundaryField(), patchI)
    {
        sigmaS_.boundaryField()[patchI] == this->sigmaS(patchI)();
    }
}


Foam::tmp<Foam::volScalarField> Foam::interpolatedSolidThermo::Cp() const
{
    tmp<volScalarField> tCp
    (
        new volScalarField
        (
            IOobject
            (
                "Cp",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimEnergy/(dimMass*dimTemperature)
        )
    );
    volScalarField& Cp = tCp();

    Cp.internalField() = interpolateXY
    (
        T_.internalField(),
        TValues_,
        cpValues_
    );

    forAll(Cp.boundaryField(), patchI)
    {
        Cp.boundaryField()[patchI] == this->Cp(patchI)();
    }

    return tCp;
}


Foam::tmp<Foam::volScalarField> Foam::interpolatedSolidThermo::Hf() const
{
    tmp<volScalarField> tHf
    (
        new volScalarField
        (
            IOobject
            (
                "Hf",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimEnergy/dimMass
        )
    );
    volScalarField& Hf = tHf();

    Hf.internalField() = interpolateXY
    (
        T_.internalField(),
        TValues_,
        HfValues_
    );

    forAll(Hf.boundaryField(), patchI)
    {
        Hf.boundaryField()[patchI] == this->Hf(patchI)();
    }

    return tHf;
}


Foam::tmp<Foam::scalarField> Foam::interpolatedSolidThermo::Cp
(
    const label patchI
) const
{
    return tmp<scalarField>
    (
        new scalarField
        (
            interpolateXY
            (
                T_.boundaryField()[patchI],
                TValues_,
                cpValues_
            )
        )
    );
}


Foam::tmp<Foam::scalarField> Foam::interpolatedSolidThermo::Hf
(
    const label patchI
) const
{
    return tmp<scalarField>
    (
        new scalarField
        (
            interpolateXY
            (
                T_.boundaryField()[patchI],
                TValues_,
                HfValues_
            )
        )
    );
}


bool Foam::interpolatedSolidThermo::read()
{
    return read(dict_);
}


bool Foam::interpolatedSolidThermo::read(const dictionary& dict)
{
    bool ok = interpolateSolid::read(dict);
    return ok;
}


bool Foam::interpolatedSolidThermo::writeData(Ostream& os) const
{
    bool ok = basicSolidThermo::writeData(os);
    ok = interpolateSolid::writeData(os);

    return ok && os.good();
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const interpolatedSolidThermo& s)
{
    s.writeData(os);
    return os;
}


// ************************************************************************* //
