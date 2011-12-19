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

#include "LocalInteraction.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class CloudType>
void Foam::LocalInteraction<CloudType>::readProps()
{
    if (!this->owner().solution().transient())
    {
        return;
    }

    IOobject propsDictHeader
    (
        "localInteractionProperties",
        this->owner().db().time().timeName(),
        "uniform"/cloud::prefix/this->owner().name(),
        this->owner().db(),
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE,
        false
    );

    if (propsDictHeader.headerOk())
    {
        const IOdictionary propsDict(propsDictHeader);
        propsDict.readIfPresent("nEscape", nEscape0_);
        propsDict.readIfPresent("massEscape", massEscape0_);
        propsDict.readIfPresent("nStick", nStick0_);
        propsDict.readIfPresent("massStick", massStick0_);
    }
}


template<class CloudType>
void Foam::LocalInteraction<CloudType>::writeProps
(
    const labelList& nEscape,
    const scalarList& massEscape,
    const labelList& nStick,
    const scalarList& massStick
) const
{
    if (!this->owner().solution().transient())
    {
        return;
    }

    if (this->owner().db().time().outputTime())
    {
        IOdictionary propsDict
        (
            IOobject
            (
                "localInteractionProperties",
                this->owner().db().time().timeName(),
                "uniform"/cloud::prefix/this->owner().name(),
                this->owner().db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        propsDict.add("nEscape", nEscape);
        propsDict.add("massEscape", massEscape);
        propsDict.add("nStick", nStick);
        propsDict.add("massStick", massStick);

        propsDict.writeObject
        (
            IOstream::ASCII,
            IOstream::currentVersion,
            this->owner().db().time().writeCompression()
        );
    }
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::LocalInteraction<CloudType>::LocalInteraction
(
    const dictionary& dict,
    CloudType& cloud
)
:
    PatchInteractionModel<CloudType>(dict, cloud, typeName),
    patchData_(cloud.mesh(), this->coeffDict()),
    nEscape0_(patchData_.size(), 0),
    massEscape0_(patchData_.size(), 0.0),
    nStick0_(patchData_.size(), 0),
    massStick0_(patchData_.size(), 0.0),
    nEscape_(patchData_.size(), 0),
    massEscape_(patchData_.size(), 0.0),
    nStick_(patchData_.size(), 0),
    massStick_(patchData_.size(), 0.0)
{
    // check that interactions are valid/specified
    forAll(patchData_, patchI)
    {
        const word& interactionTypeName =
            patchData_[patchI].interactionTypeName();
        const typename PatchInteractionModel<CloudType>::interactionType& it =
            this->wordToInteractionType(interactionTypeName);

        if (it == PatchInteractionModel<CloudType>::itOther)
        {
            const word& patchName = patchData_[patchI].patchName();
            FatalErrorIn("LocalInteraction(const dictionary&, CloudType&)")
                << "Unknown patch interaction type "
                << interactionTypeName << " for patch " << patchName
                << ". Valid selections are:"
                << this->PatchInteractionModel<CloudType>::interactionTypeNames_
                << nl << exit(FatalError);
        }
    }

    readProps();
}


template<class CloudType>
Foam::LocalInteraction<CloudType>::LocalInteraction
(
    const LocalInteraction<CloudType>& pim
)
:
    PatchInteractionModel<CloudType>(pim),
    patchData_(pim.patchData_),
    nEscape0_(pim.nEscape0_),
    massEscape0_(pim.massEscape0_),
    nStick0_(pim.nStick0_),
    massStick0_(pim.massStick0_),
    nEscape_(pim.nEscape_),
    massEscape_(pim.massEscape_),
    nStick_(pim.nStick_),
    massStick_(pim.massStick_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::LocalInteraction<CloudType>::~LocalInteraction()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::LocalInteraction<CloudType>::correct
(
    typename CloudType::parcelType& p,
    const polyPatch& pp,
    bool& keepParticle,
    const scalar trackFraction,
    const tetIndices& tetIs
)
{
    vector& U = p.U();

    bool& active = p.active();

    label patchI = patchData_.applyToPatch(pp.index());

    if (patchI >= 0)
    {
        typename PatchInteractionModel<CloudType>::interactionType it =
            this->wordToInteractionType
            (
                patchData_[patchI].interactionTypeName()
            );

        switch (it)
        {
            case PatchInteractionModel<CloudType>::itEscape:
            {
                keepParticle = false;
                active = false;
                U = vector::zero;
                nEscape_[patchI]++;
                massEscape_[patchI] += p.mass()*p.nParticle();
                break;
            }
            case PatchInteractionModel<CloudType>::itStick:
            {
                keepParticle = true;
                active = false;
                U = vector::zero;
                nStick_[patchI]++;
                massStick_[patchI] += p.mass()*p.nParticle();
                break;
            }
            case PatchInteractionModel<CloudType>::itRebound:
            {
                keepParticle = true;
                active = true;

                vector nw;
                vector Up;

                this->patchData(p, pp, trackFraction, tetIs, nw, Up);

                // Calculate motion relative to patch velocity
                U -= Up;

                scalar Un = U & nw;
                vector Ut = U - Un*nw;

                if (Un > 0)
                {
                    U -= (1.0 + patchData_[patchI].e())*Un*nw;
                }

                U -= patchData_[patchI].mu()*Ut;

                // Return velocity to global space
                U += Up;

                break;
            }
            default:
            {
                FatalErrorIn
                (
                    "bool LocalInteraction<CloudType>::correct"
                    "("
                        "typename CloudType::parcelType&, "
                        "const polyPatch&, "
                        "bool&, "
                        "scalar&, "
                        "const tetIndices&"
                    ") const"
                )   << "Unknown interaction type "
                    << patchData_[patchI].interactionTypeName()
                    << "(" << it << ") for patch "
                    << patchData_[patchI].patchName()
                    << ". Valid selections are:" << this->interactionTypeNames_
                    << endl << abort(FatalError);
            }
        }

        return true;
    }

    return false;
}


template<class CloudType>
void Foam::LocalInteraction<CloudType>::info(Ostream& os) const
{
    labelList npe(nEscape_);
    Pstream::listCombineGather(npe, plusEqOp<label>());
    npe = npe + nEscape0_;

    scalarList mpe(massEscape_);
    Pstream::listCombineGather(mpe, plusEqOp<scalar>());
    mpe = mpe + massEscape0_;

    labelList nps(nStick_);
    Pstream::listCombineGather(nps, plusEqOp<label>());
    nps = nps + nStick0_;

    scalarList mps(massStick_);
    Pstream::listCombineGather(mps, plusEqOp<scalar>());
    mps = mps + massStick0_;


    forAll(patchData_, i)
    {
        os  << "    Parcel fate (number, mass)      : patch "
            <<  patchData_[i].patchName() << nl
            << "      - escape                      = " << npe[i]
            << ", " << mpe[i] << nl
            << "      - stick                       = " << nps[i]
            << ", " << mps[i] << nl;
    }

    writeProps(npe, mpe, nps, mps);
}


// ************************************************************************* //
