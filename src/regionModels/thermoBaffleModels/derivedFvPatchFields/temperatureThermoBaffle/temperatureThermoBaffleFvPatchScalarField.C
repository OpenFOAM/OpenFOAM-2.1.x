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

#include "temperatureThermoBaffleFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

temperatureThermoBaffleFvPatchScalarField::
temperatureThermoBaffleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    turbulentTemperatureCoupledBaffleMixedFvPatchScalarField(p, iF),
    owner_(false),
    baffle_(),
    solidThermoType_("undefined"),
    dict_(dictionary::null)
{}


temperatureThermoBaffleFvPatchScalarField::
temperatureThermoBaffleFvPatchScalarField
(
    const temperatureThermoBaffleFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    turbulentTemperatureCoupledBaffleMixedFvPatchScalarField
    (
        ptf,
        p,
        iF,
        mapper
    ),
    owner_(ptf.owner_),
    baffle_(ptf.baffle_),
    solidThermoType_(ptf.solidThermoType_),
    dict_(ptf.dict_)
{}


temperatureThermoBaffleFvPatchScalarField::
temperatureThermoBaffleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    turbulentTemperatureCoupledBaffleMixedFvPatchScalarField(p, iF, dict),
    owner_(false),
    baffle_(),
    solidThermoType_(),
    dict_(dict)
{
    if (!isA<mappedPatchBase>(patch().patch()))
    {
        FatalErrorIn
        (
            "temperatureThermoBaffleFvPatchScalarField::"
            "temperatureThermoBaffleFvPatchScalarField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<scalar, volMesh>& iF,\n"
            "    const dictionary& dict\n"
            ")\n"
        )   << "\n    patch type '" << patch().type()
            << "' not type '" << mappedPatchBase::typeName << "'"
            << "\n    for patch " << patch().name()
            << " of field " << dimensionedInternalField().name()
            << " in file " << dimensionedInternalField().objectPath()
            << exit(FatalError);
    }

    const fvMesh& thisMesh = patch().boundaryMesh().mesh();

    typedef regionModels::thermoBaffleModels::thermoBaffleModel baffle;

    if
    (
        thisMesh.name() == polyMesh::defaultRegion
     && !thisMesh.foundObject<baffle>("thermoBaffle")
     && !owner_
    )
    {
        Info << "Creating thermal baffle..." << endl;
        baffle_.reset(baffle::New(thisMesh, dict).ptr());
        owner_ = true;
        dict.lookup("thermoType") >> solidThermoType_;
    }

}


temperatureThermoBaffleFvPatchScalarField::
temperatureThermoBaffleFvPatchScalarField
(
    const temperatureThermoBaffleFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    turbulentTemperatureCoupledBaffleMixedFvPatchScalarField(ptf, iF),
    owner_(ptf.owner_),
    baffle_(ptf.baffle_),
    solidThermoType_(ptf.solidThermoType_),
    dict_(ptf.dict_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void temperatureThermoBaffleFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);
}


void temperatureThermoBaffleFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);
}


void temperatureThermoBaffleFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const fvMesh& thisMesh = patch().boundaryMesh().mesh();

    if
    (
        thisMesh.name() == polyMesh::defaultRegion
     && owner_
    )
    {
        baffle_->evolve();
    }

    turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::updateCoeffs();
}


void temperatureThermoBaffleFvPatchScalarField::write(Ostream& os) const
{
    turbulentTemperatureCoupledBaffleMixedFvPatchScalarField::write(os);

    const fvMesh& thisMesh = patch().boundaryMesh().mesh();

    if (thisMesh.name() == polyMesh::defaultRegion && owner_)
    {
        word thermoModel = dict_.lookup("thermoBaffleModel");

        os.writeKeyword("thermoBaffleModel")
            << thermoModel
            << token::END_STATEMENT << nl;

        word regionName = dict_.lookup("regionName");
        os.writeKeyword("regionName") << regionName
            << token::END_STATEMENT << nl;

        bool infoOutput = readBool(dict_.lookup("infoOutput"));
        os.writeKeyword("infoOutput") << infoOutput
            << token::END_STATEMENT << nl;

        bool active = readBool(dict_.lookup("active"));
        os.writeKeyword("active") <<  active
            << token::END_STATEMENT << nl;

        os.writeKeyword(word(thermoModel + "Coeffs"));
        os << dict_.subDict(thermoModel + "Coeffs") << nl;

        os.writeKeyword("thermoType") << solidThermoType_
            << token::END_STATEMENT << nl;

        os.writeKeyword(word(solidThermoType_ + "Coeffs"));

        os << dict_.subDict(solidThermoType_ + "Coeffs") << nl;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    temperatureThermoBaffleFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
