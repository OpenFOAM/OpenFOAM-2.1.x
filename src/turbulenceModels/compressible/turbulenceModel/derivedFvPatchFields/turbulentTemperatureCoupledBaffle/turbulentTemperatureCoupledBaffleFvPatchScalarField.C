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

#include "turbulentTemperatureCoupledBaffleFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"
#include "regionProperties.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::turbulentTemperatureCoupledBaffleFvPatchScalarField::interfaceOwner
(
    const polyMesh& nbrRegion,
    const polyPatch& nbrPatch
) const
{
    const fvMesh& myRegion = patch().boundaryMesh().mesh();

    if (nbrRegion.name() == myRegion.name())
    {
        return patch().index() < nbrPatch.index();
    }
    else
    {
        const regionProperties& props =
            myRegion.objectRegistry::parent().lookupObject<regionProperties>
            (
                "regionProperties"
            );

        label myIndex = findIndex(props.fluidRegionNames(), myRegion.name());
        if (myIndex == -1)
        {
            label i = findIndex(props.solidRegionNames(), myRegion.name());

            if (i == -1)
            {
                FatalErrorIn
                (
                    "turbulentTemperatureCoupledBaffleFvPatchScalarField"
                    "::interfaceOwner(const polyMesh&"
                    ", const polyPatch&)const"
                )   << "Cannot find region " << myRegion.name()
                    << " neither in fluids " << props.fluidRegionNames()
                    << " nor in solids " << props.solidRegionNames()
                    << exit(FatalError);
            }
            myIndex = props.fluidRegionNames().size() + i;
        }
        label nbrIndex = findIndex
        (
            props.fluidRegionNames(),
            nbrRegion.name()
        );
        if (nbrIndex == -1)
        {
            label i = findIndex(props.solidRegionNames(), nbrRegion.name());

            if (i == -1)
            {
                FatalErrorIn
                (
                    "coupleManager::interfaceOwner"
                    "(const polyMesh&, const polyPatch&) const"
                )   << "Cannot find region " << nbrRegion.name()
                    << " neither in fluids " << props.fluidRegionNames()
                    << " nor in solids " << props.solidRegionNames()
                    << exit(FatalError);
            }
            nbrIndex = props.fluidRegionNames().size() + i;
        }

        return myIndex < nbrIndex;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbulentTemperatureCoupledBaffleFvPatchScalarField::
turbulentTemperatureCoupledBaffleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), "undefined", "undefined-K"),
    neighbourFieldName_("undefined-neighbourFieldName")
{}


Foam::turbulentTemperatureCoupledBaffleFvPatchScalarField::
turbulentTemperatureCoupledBaffleFvPatchScalarField
(
    const turbulentTemperatureCoupledBaffleFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    temperatureCoupledBase(patch(), ptf.KMethod(), ptf.KName()),
    neighbourFieldName_(ptf.neighbourFieldName_)
{}


Foam::turbulentTemperatureCoupledBaffleFvPatchScalarField::
turbulentTemperatureCoupledBaffleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    temperatureCoupledBase(patch(), dict),
    neighbourFieldName_(dict.lookup("neighbourFieldName"))
{
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalErrorIn
        (
            "turbulentTemperatureCoupledBaffleFvPatchScalarField::"
            "turbulentTemperatureCoupledBaffleFvPatchScalarField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<scalar, volMesh>& iF,\n"
            "    const dictionary& dict\n"
            ")\n"
        )   << "\n    patch type '" << p.type()
            << "' not type '" << mappedPatchBase::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << dimensionedInternalField().name()
            << " in file " << dimensionedInternalField().objectPath()
            << exit(FatalError);
    }
}


Foam::turbulentTemperatureCoupledBaffleFvPatchScalarField::
turbulentTemperatureCoupledBaffleFvPatchScalarField
(
    const turbulentTemperatureCoupledBaffleFvPatchScalarField& wtcsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(wtcsf, iF),
    temperatureCoupledBase(patch(), wtcsf.KMethod(), wtcsf.KName()),
    neighbourFieldName_(wtcsf.neighbourFieldName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::turbulentTemperatureCoupledBaffleFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    // Get the coupling information from the mappedPatchBase
    const mappedPatchBase& mpp = refCast<const mappedPatchBase>
    (
        this->patch().patch()
    );
    const polyMesh& nbrMesh = mpp.sampleMesh();
    const fvPatch& nbrPatch = refCast<const fvMesh>
    (
        nbrMesh
    ).boundary()[mpp.samplePolyPatch().index()];

    // Force recalculation of mapping and schedule
    const mapDistribute& distMap = mpp.map();
    (void)distMap.schedule();

    tmp<scalarField> intFld = patchInternalField();

    if (interfaceOwner(nbrMesh, nbrPatch.patch()))
    {
        // Note: other side information could be cached - it only needs
        // to be updated the first time round the iteration (i.e. when
        // switching regions) but unfortunately we don't have this information.


        // Calculate the temperature by harmonic averaging
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        const turbulentTemperatureCoupledBaffleFvPatchScalarField& nbrField =
        refCast<const turbulentTemperatureCoupledBaffleFvPatchScalarField>
        (
            nbrPatch.lookupPatchField<volScalarField, scalar>
            (
                neighbourFieldName_
            )
        );

        // Swap to obtain full local values of neighbour internal field
        scalarField nbrIntFld(nbrField.patchInternalField());
        distMap.distribute(nbrIntFld);

        // Swap to obtain full local values of neighbour K*delta
        scalarField nbrKDelta(nbrField.K(nbrField)*nbrPatch.deltaCoeffs());
        distMap.distribute(nbrKDelta);

        tmp<scalarField> myKDelta = K(*this)*patch().deltaCoeffs();

        // Calculate common wall temperature. Reuse *this to store common value.
        scalarField Twall
        (
            (myKDelta()*intFld() + nbrKDelta*nbrIntFld)
          / (myKDelta() + nbrKDelta)
        );
        // Assign to me
        fvPatchScalarField::operator=(Twall);
        // Distribute back and assign to neighbour
        distMap.reverseDistribute(nbrField.size(), Twall);
        const_cast<turbulentTemperatureCoupledBaffleFvPatchScalarField&>
        (
            nbrField
        ).fvPatchScalarField::operator=(Twall);
    }

    if (debug)
    {
        //tmp<scalarField> normalGradient =
        //    (*this-intFld())
        //  * patch().deltaCoeffs();

        scalar Q = gSum(K(*this)*patch().magSf()*snGrad());

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->dimensionedInternalField().name() << " <- "
            << nbrMesh.name() << ':'
            << nbrPatch.name() << ':'
            << this->dimensionedInternalField().name() << " :"
            << " heat[W]:" << Q
            << " walltemperature "
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << endl;
    }

    // Restore tag
    UPstream::msgType() = oldTag;

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::turbulentTemperatureCoupledBaffleFvPatchScalarField::write
(
    Ostream& os
) const
{
    fixedValueFvPatchScalarField::write(os);
    os.writeKeyword("neighbourFieldName")<< neighbourFieldName_
        << token::END_STATEMENT << nl;
    temperatureCoupledBase::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

makePatchTypeField
(
    fvPatchScalarField,
    turbulentTemperatureCoupledBaffleFvPatchScalarField
);

} // End namespace Foam

// ************************************************************************* //
