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

#include "porousZone.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "geometricOneField.H"

// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

// adjust negative resistance values to be multiplier of max value
void Foam::porousZone::adjustNegativeResistance(dimensionedVector& resist)
{
    scalar maxCmpt = max(0, cmptMax(resist.value()));

    if (maxCmpt < 0)
    {
        FatalErrorIn
        (
            "Foam::porousZone::porousZone::adjustNegativeResistance"
            "(dimensionedVector&)"
        )   << "negative resistances! " << resist
            << exit(FatalError);
    }
    else
    {
        vector& val = resist.value();
        for (label cmpt=0; cmpt < vector::nComponents; ++cmpt)
        {
            if (val[cmpt] < 0)
            {
                val[cmpt] *= -maxCmpt;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porousZone::porousZone
(
    const keyType& key,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    key_(key),
    mesh_(mesh),
    dict_(dict),
    cellZoneIds_(mesh_.cellZones().findIndices(key_)),
    coordSys_(dict, mesh),
    porosity_(1),
    intensity_(0),
    mixingLength_(0),
    C0_(0),
    C1_(0),
    D_("D", dimensionSet(0, -2, 0, 0, 0), tensor::zero),
    F_("F", dimensionSet(0, -1, 0, 0, 0), tensor::zero)
{
    Info<< "Creating porous zone: " << key_ << endl;

    bool foundZone = !cellZoneIds_.empty();
    reduce(foundZone, orOp<bool>());

    if (!foundZone && Pstream::master())
    {
        FatalErrorIn
        (
            "Foam::porousZone::porousZone"
            "(const keyType&, const fvMesh&, const dictionary&)"
        )   << "cannot find porous cellZone " << key_
            << exit(FatalError);
    }


    // porosity
    if
    (
        dict_.readIfPresent("porosity", porosity_)
     && (porosity_ <= 0.0 || porosity_ > 1.0)
    )
    {
        FatalIOErrorIn
        (
            "Foam::porousZone::porousZone"
            "(const keyType&, const fvMesh&, const dictionary&)",
            dict_
        )
            << "out-of-range porosity value " << porosity_
            << exit(FatalIOError);
    }

    // turbulent intensity
    if
    (
        dict_.readIfPresent("intensity", intensity_)
     && (intensity_ <= 0.0 || intensity_ > 1.0)
    )
    {
        FatalIOErrorIn
        (
            "Foam::porousZone::porousZone"
            "(const keyType&, const fvMesh&, const dictionary&)",
            dict_
        )
            << "out-of-range turbulent intensity value " << intensity_
            << exit(FatalIOError);
    }

    // turbulent length scale
    if
    (
        dict_.readIfPresent("mixingLength", mixingLength_)
     && (mixingLength_ <= 0.0)
    )
    {
        FatalIOErrorIn
        (
            "Foam::porousZone::porousZone"
            "(const keyType&, const fvMesh&, const dictionary&)",
            dict_
        )
            << "out-of-range turbulent length scale " << mixingLength_
            << exit(FatalIOError);
    }


    // powerLaw coefficients
    if (const dictionary* dictPtr = dict_.subDictPtr("powerLaw"))
    {
        dictPtr->readIfPresent("C0", C0_);
        dictPtr->readIfPresent("C1", C1_);
    }

    // Darcy-Forchheimer coefficients
    if (const dictionary* dictPtr = dict_.subDictPtr("Darcy"))
    {
        // local-to-global transformation tensor
        const tensor& E = coordSys_.R();

        dimensionedVector d(vector::zero);
        if (dictPtr->readIfPresent("d", d))
        {
            if (D_.dimensions() != d.dimensions())
            {
                FatalIOErrorIn
                (
                    "Foam::porousZone::porousZone"
                    "(const keyType&, const fvMesh&, const dictionary&)",
                    dict_
                )   << "incorrect dimensions for d: " << d.dimensions()
                    << " should be " << D_.dimensions()
                    << exit(FatalIOError);
            }

            adjustNegativeResistance(d);

            D_.value().xx() = d.value().x();
            D_.value().yy() = d.value().y();
            D_.value().zz() = d.value().z();
            D_.value() = (E & D_ & E.T()).value();
        }

        dimensionedVector f(vector::zero);
        if (dictPtr->readIfPresent("f", f))
        {
            if (F_.dimensions() != f.dimensions())
            {
                FatalIOErrorIn
                (
                    "Foam::porousZone::porousZone"
                    "(const keyType&, const fvMesh&, const dictionary&)",
                    dict_
                )   << "incorrect dimensions for f: " << f.dimensions()
                    << " should be " << F_.dimensions()
                    << exit(FatalIOError);
            }

            adjustNegativeResistance(f);

            // leading 0.5 is from 1/2 * rho
            F_.value().xx() = 0.5*f.value().x();
            F_.value().yy() = 0.5*f.value().y();
            F_.value().zz() = 0.5*f.value().z();
            F_.value() = (E & F_ & E.T()).value();
        }
    }

    // it is an error not to define anything
    if
    (
        C0_ <= VSMALL
     && magSqr(D_.value()) <= VSMALL
     && magSqr(F_.value()) <= VSMALL
    )
    {
        FatalIOErrorIn
        (
            "Foam::porousZone::porousZone"
            "(const keyType&, const fvMesh&, const dictionary&)",
            dict_
        )   << "neither powerLaw (C0/C1) "
               "nor Darcy-Forchheimer law (d/f) specified"
            << exit(FatalIOError);
    }

    // feedback for the user
    if (dict.lookupOrDefault("printCoeffs", false))
    {
        writeDict(Info, false);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::porousZone::addResistance(fvVectorMatrix& UEqn) const
{
    if (cellZoneIds_.empty())
    {
        return;
    }

    bool compressible = false;
    if (UEqn.dimensions() == dimensionSet(1, 1, -2, 0, 0))
    {
        compressible = true;
    }

    const scalarField& V = mesh_.V();
    scalarField& Udiag = UEqn.diag();
    vectorField& Usource = UEqn.source();
    const vectorField& U = UEqn.psi();

    if (C0_ > VSMALL)
    {
        if (compressible)
        {
            addPowerLawResistance
            (
                Udiag,
                V,
                mesh_.lookupObject<volScalarField>("rho"),
                U
            );
        }
        else
        {
            addPowerLawResistance
            (
                Udiag,
                V,
                geometricOneField(),
                U
            );
        }
    }

    const tensor& D = D_.value();
    const tensor& F = F_.value();

    if (magSqr(D) > VSMALL || magSqr(F) > VSMALL)
    {
        if (compressible)
        {
            addViscousInertialResistance
            (
                Udiag,
                Usource,
                V,
                mesh_.lookupObject<volScalarField>("rho"),
                mesh_.lookupObject<volScalarField>("mu"),
                U
            );
        }
        else
        {
            addViscousInertialResistance
            (
                Udiag,
                Usource,
                V,
                geometricOneField(),
                mesh_.lookupObject<volScalarField>("nu"),
                U
            );
        }
    }
}


void Foam::porousZone::addResistance
(
    fvVectorMatrix& UEqn,
    const volScalarField& rho,
    const volScalarField& mu
) const
{
    if (cellZoneIds_.empty())
    {
        return;
    }

    const scalarField& V = mesh_.V();
    scalarField& Udiag = UEqn.diag();
    vectorField& Usource = UEqn.source();
    const vectorField& U = UEqn.psi();

    if (C0_ > VSMALL)
    {
        addPowerLawResistance
        (
            Udiag,
            V,
            rho,
            U
        );
    }

    const tensor& D = D_.value();
    const tensor& F = F_.value();

    if (magSqr(D) > VSMALL || magSqr(F) > VSMALL)
    {
        addViscousInertialResistance
        (
            Udiag,
            Usource,
            V,
            rho,
            mu,
            U
        );
    }
}

void Foam::porousZone::addResistance
(
    const fvVectorMatrix& UEqn,
    volTensorField& AU,
    bool correctAUprocBC
) const
{
    if (cellZoneIds_.empty())
    {
        return;
    }

    bool compressible = false;
    if (UEqn.dimensions() == dimensionSet(1, 1, -2, 0, 0))
    {
        compressible = true;
    }

    const vectorField& U = UEqn.psi();

    if (C0_ > VSMALL)
    {
        if (compressible)
        {
            addPowerLawResistance
            (
                AU,
                mesh_.lookupObject<volScalarField>("rho"),
                U
            );
        }
        else
        {
            addPowerLawResistance
            (
                AU,
                geometricOneField(),
                U
            );
        }
    }

    const tensor& D = D_.value();
    const tensor& F = F_.value();

    if (magSqr(D) > VSMALL || magSqr(F) > VSMALL)
    {
        if (compressible)
        {
            addViscousInertialResistance
            (
                AU,
                mesh_.lookupObject<volScalarField>("rho"),
                mesh_.lookupObject<volScalarField>("mu"),
                U
            );
        }
        else
        {
            addViscousInertialResistance
            (
                AU,
                geometricOneField(),
                mesh_.lookupObject<volScalarField>("nu"),
                U
            );
        }
    }

    if (correctAUprocBC)
    {
        // Correct the boundary conditions of the tensorial diagonal to ensure
        // processor boundaries are correctly handled when AU^-1 is interpolated
        // for the pressure equation.
        AU.correctBoundaryConditions();
    }
}


void Foam::porousZone::writeDict(Ostream& os, bool subDict) const
{
    if (subDict)
    {
        os  << indent << token::BEGIN_BLOCK << incrIndent << nl;
        os.writeKeyword("name")
            << zoneName() << token::END_STATEMENT << nl;
    }
    else
    {
        os  << indent << zoneName() << nl
            << indent << token::BEGIN_BLOCK << incrIndent << nl;
    }

    if (dict_.found("note"))
    {
        os.writeKeyword("note")
            << string(dict_.lookup("note")) << token::END_STATEMENT << nl;
    }

    coordSys_.writeDict(os, true);

    if (dict_.found("porosity"))
    {
        os.writeKeyword("porosity")
            << porosity() << token::END_STATEMENT << nl;
    }

    if (dict_.found("intensity"))
    {
        os.writeKeyword("intensity")
            << intensity() << token::END_STATEMENT << nl;
    }

    if (dict_.found("mixingLength"))
    {
        os.writeKeyword("mixingLength")
            << mixingLength() << token::END_STATEMENT << nl;
    }

    // powerLaw coefficients
    if (const dictionary* dictPtr = dict_.subDictPtr("powerLaw"))
    {
        os  << indent << "powerLaw";
        dictPtr->write(os);
    }

    // Darcy-Forchheimer coefficients
    if (const dictionary* dictPtr = dict_.subDictPtr("Darcy"))
    {
        os  << indent << "Darcy";
        dictPtr->write(os);
    }

    // thermalModel
    if (const dictionary* dictPtr = dict_.subDictPtr("thermalModel"))
    {
        os  << indent << "thermalModel";
        dictPtr->write(os);
    }

    os  << decrIndent << indent << token::END_BLOCK << endl;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const porousZone& pz)
{
    pz.writeDict(os);
    return os;
}

// ************************************************************************* //
