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

#include "temperatureCoupledBase.H"
#include "volFields.H"
#include "basicSolidThermo.H"
#include "turbulenceModel.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* Foam::NamedEnum
    <
        Foam::temperatureCoupledBase::KMethodType,
        4
    >::names[] =
    {
        "basicThermo",
        "solidThermo",
        "directionalSolidThermo",
        "lookup"
    };
}


const Foam::NamedEnum<Foam::temperatureCoupledBase::KMethodType, 4>
    Foam::temperatureCoupledBase::KMethodTypeNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::temperatureCoupledBase::temperatureCoupledBase
(
    const fvPatch& patch,
    const word& calculationType,
    const word& KName
)
:
    patch_(patch),
    method_(KMethodTypeNames_[calculationType]),
    KName_(KName)
{}


Foam::temperatureCoupledBase::temperatureCoupledBase
(
    const fvPatch& patch,
    const dictionary& dict
)
:
    patch_(patch),
    method_(KMethodTypeNames_.read(dict.lookup("K"))),
    KName_(dict.lookup("KName"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::temperatureCoupledBase::K
(
    const scalarField& Tp
) const
{
    const fvMesh& mesh = patch_.boundaryMesh().mesh();

    switch (method_)
    {
        case BASICTHERMO:
        {
            const compressible::turbulenceModel& model =
                mesh.lookupObject<compressible::turbulenceModel>
                (
                    "turbulenceModel"
                );

            return
                model.alphaEff()().boundaryField()[patch_.index()]
               *model.thermo().Cp(Tp, patch_.index());
        }
        break;

        case SOLIDTHERMO:
        {
            const basicSolidThermo& thermo =
                mesh.lookupObject<basicSolidThermo>
                (
                    "solidThermophysicalProperties"
                );
            return thermo.K(patch_.index());
        }
        break;

        case DIRECTIONALSOLIDTHERMO:
        {
            const vectorField n(patch_.nf());

            const basicSolidThermo& thermo =
                mesh.lookupObject<basicSolidThermo>
                (
                    "solidThermophysicalProperties"
                );
            return n & thermo.directionalK(patch_.index()) & n;
        }
        break;

        case LOOKUP:
        {
            if (mesh.objectRegistry::foundObject<volScalarField>(KName_))
            {
                return patch_.lookupPatchField<volScalarField, scalar>(KName_);
            }
            else if
            (
                mesh.objectRegistry::foundObject<volSymmTensorField>(KName_)
            )
            {
                const symmTensorField& KWall =
                    patch_.lookupPatchField<volSymmTensorField, scalar>(KName_);

                const vectorField n(patch_.nf());

                return n & KWall & n;
            }
            else
            {
                FatalErrorIn("temperatureCoupledBase::K() const")
                    << "Did not find field " << KName_
                    << " on mesh " << mesh.name() << " patch " << patch_.name()
                    << endl
                    << "Please set 'K' to one of " << KMethodTypeNames_.toc()
                    << " and 'KName' to the name of the volScalar"
                    << " or volSymmTensor field (if K=lookup)"
                    << exit(FatalError);

                return scalarField(0);
            }
        }

        default:
        {
            FatalErrorIn("temperatureCoupledBase::K() const")
                << "Unimplemented method " << method_ << endl
                << "Please set 'K' to one of " << KMethodTypeNames_.toc()
                << " and 'KName' to the name of the volScalar"
                << " or volSymmTensor field (if K=lookup)"
                << exit(FatalError);
        }
        break;
    }
    return scalarField(0);
}


void Foam::temperatureCoupledBase::write(Ostream& os) const
{
    os.writeKeyword("K") << KMethodTypeNames_[method_]
        << token::END_STATEMENT << nl;
    os.writeKeyword("KName") << KName_ << token::END_STATEMENT << nl;
}


// ************************************************************************* //
