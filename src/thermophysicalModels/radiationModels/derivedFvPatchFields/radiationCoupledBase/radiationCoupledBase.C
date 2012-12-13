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

#include "radiationCoupledBase.H"
#include "volFields.H"
#include "basicSolidThermo.H"

#include "mappedPatchBase.H"
#include "fvPatchFieldMapper.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* Foam::NamedEnum
    <
        Foam::radiationCoupledBase::emissivityMethodType,
        2
    >::names[] =
    {
        "solidThermo",
        "lookup"
    };
}


const Foam::NamedEnum<Foam::radiationCoupledBase::emissivityMethodType, 2>
    Foam::radiationCoupledBase::emissivityMethodTypeNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationCoupledBase::radiationCoupledBase
(
    const fvPatch& patch,
    const word& calculationType,
    const scalarField& emissivity
)
:
    patch_(patch),
    method_(emissivityMethodTypeNames_[calculationType]),
    emissivity_(emissivity)
{}


Foam::radiationCoupledBase::radiationCoupledBase
(
    const fvPatch& patch,
    const dictionary& dict
)
:
    patch_(patch),
    method_(emissivityMethodTypeNames_.read(dict.lookup("emissivityMode")))
{
    switch (method_)
    {
        case SOLIDTHERMO:
        {
            if (!isA<mappedPatchBase>(patch_.patch()))
            {
                FatalIOErrorIn
                (
                    "radiationCoupledBase::radiationCoupledBase\n"
                    "(\n"
                    "    const fvPatch& p,\n"
                    "    const dictionary& dict\n"
                    ")\n",
                    dict
                )   << "\n    patch type '" << patch_.type()
                    << "' not type '" << mappedPatchBase::typeName << "'"
                    << "\n    for patch " << patch_.name()
                    << exit(FatalIOError);
            }

            emissivity_ = scalarField(patch_.size(), 0.0);
        }
        break;

        case LOOKUP:
        {
            if(!dict.found("emissivity"))
            {
                FatalIOErrorIn
                (
                    "radiationCoupledBase::radiationCoupledBase\n"
                    "(\n"
                    "    const fvPatch& p,\n"
                    "    const dictionary& dict\n"
                    ")\n",
                    dict
                )   << "\n    emissivity key does not exist for patch "
                    << patch_.name()
                    << exit(FatalIOError);
            }
            else
            {
                emissivity_ = scalarField("emissivity", dict, patch_.size());
            }
        }
        break;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalarField Foam::radiationCoupledBase::emissivity() const
{
    switch (method_)
    {
        case SOLIDTHERMO:
        {
            // Get the coupling information from the mappedPatchBase
            const mappedPatchBase& mpp =
                refCast<const mappedPatchBase>
                (
                    patch_.patch()
                );

            const polyMesh& nbrMesh = mpp.sampleMesh();

            const fvPatch& nbrPatch = refCast<const fvMesh>
            (
                nbrMesh
            ).boundary()[mpp.samplePolyPatch().index()];

            if (nbrMesh.foundObject<volScalarField>("emissivity"))
            {
                tmp<scalarField> temissivity
                (
                    new scalarField
                    (
                        nbrPatch.lookupPatchField<volScalarField, scalar>
                        (
                            "emissivity"
                        )
                    )
                );

                scalarField emissivity(temissivity);
                // Use direct map mapping to exchange data
                mpp.distribute(emissivity);
                //Pout << emissivity << endl;
                return emissivity;
            }
            else
            {
                return scalarField(0);
            }

        }
        break;

        case LOOKUP:
        {
            // return local value
            return emissivity_;
        }

        default:
        {
            FatalErrorIn
            (
                "radiationCoupledBase::emissivity(const scalarField&)"
            )
                << "Unimplemented method " << method_ << endl
                << "Please set 'emissivity' to one of "
                << emissivityMethodTypeNames_.toc()
                << " and 'emissivityName' to the name of the volScalar"
                << exit(FatalError);
        }
        break;
    }
    return scalarField(0);
}


void Foam::radiationCoupledBase::write(Ostream& os) const
{
    os.writeKeyword("emissivityMode") << emissivityMethodTypeNames_[method_]
        << token::END_STATEMENT << nl;
    emissivity_.writeEntry("emissivity", os);
}


// ************************************************************************* //
