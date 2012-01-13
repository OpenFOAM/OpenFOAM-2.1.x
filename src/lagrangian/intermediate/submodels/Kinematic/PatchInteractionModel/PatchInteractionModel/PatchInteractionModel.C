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

#include "PatchInteractionModel.H"
#include "fvMesh.H"
#include "Time.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class CloudType>
Foam::wordList Foam::PatchInteractionModel<CloudType>::interactionTypeNames_
(
    IStringStream
    (
        "(rebound stick escape)"
    )()
);

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::word Foam::PatchInteractionModel<CloudType>::interactionTypeToWord
(
    const interactionType& itEnum
)
{
    word it = "other";

    switch (itEnum)
    {
        case itRebound:
        {
            it = "rebound";
            break;
        }
        case itStick:
        {
            it = "stick";
            break;
        }
        case itEscape:
        {
            it = "escape";
            break;
        }
        default:
        {
        }
    }

    return it;
}


template<class CloudType>
typename Foam::PatchInteractionModel<CloudType>::interactionType
Foam::PatchInteractionModel<CloudType>::wordToInteractionType
(
    const word& itWord
)
{
    if (itWord == "rebound")
    {
        return itRebound;
    }
    else if (itWord == "stick")
    {
        return itStick;
    }
    else if (itWord == "escape")
    {
        return itEscape;
    }
    else
    {
        return itOther;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PatchInteractionModel<CloudType>::PatchInteractionModel
(
    CloudType& owner
)
:
    SubModelBase<CloudType>(owner),
    UName_("unknown_UName")
{}


template<class CloudType>
Foam::PatchInteractionModel<CloudType>::PatchInteractionModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
:
    SubModelBase<CloudType>(owner, dict, typeName, type),
    UName_(this->coeffDict().lookupOrDefault("UName", word("U")))
{}


template<class CloudType>
Foam::PatchInteractionModel<CloudType>::PatchInteractionModel
(
    const PatchInteractionModel<CloudType>& pim
)
:
    SubModelBase<CloudType>(pim),
    UName_(pim.UName_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PatchInteractionModel<CloudType>::~PatchInteractionModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
const Foam::word& Foam::PatchInteractionModel<CloudType>::UName() const
{
    return UName_;
}


template<class CloudType>
bool Foam::PatchInteractionModel<CloudType>::correct
(
    typename CloudType::parcelType&,
    const polyPatch&,
    bool&,
    const scalar,
    const tetIndices&
)
{
    notImplemented
    (
        "bool Foam::PatchInteractionModel<CloudType>::correct"
        "("
            "typename CloudType::parcelType&, "
            "const polyPatch&, "
            "bool&, "
            "const scalar, "
            "const tetIndices& "
        ") const"
    );
    return false;
}


template<class CloudType>
void Foam::PatchInteractionModel<CloudType>::patchData
(
    typename CloudType::parcelType& p,
    const polyPatch& pp,
    const scalar trackFraction,
    const tetIndices& tetIs,
    vector& nw,
    vector& Up
) const
{
    const fvMesh& mesh = this->owner().mesh();

    const volVectorField& Ufield =
        mesh.objectRegistry::lookupObject<volVectorField>(UName_);

    label patchI = pp.index();
    label patchFaceI = pp.whichFace(p.face());

    vector n = tetIs.faceTri(mesh).normal();
    n /= mag(n);

    vector U = Ufield.boundaryField()[patchI][patchFaceI];

    // Unless the face is rotating, the required normal is n;
    nw = n;

    if (!mesh.moving())
    {
        // Only wall patches may have a non-zero wall velocity from
        // the velocity field when the mesh is not moving.

        if (isA<wallPolyPatch>(pp))
        {
            Up = U;
        }
        else
        {
            Up = vector::zero;
        }
    }
    else
    {
        vector U00 = Ufield.oldTime().boundaryField()[patchI][patchFaceI];

        vector n00 = tetIs.oldFaceTri(mesh).normal();

        // Difference in normal over timestep
        vector dn = vector::zero;

        if (mag(n00) > SMALL)
        {
            // If the old normal is zero (for example in layer
            // addition) then use the current normal, meaning that the
            // motion can only be translational, and dn remains zero,
            // otherwise, calculate dn:

            n00 /= mag(n00);

            dn = n - n00;
        }

        // Total fraction thought the timestep of the motion,
        // including stepFraction before the current tracking step
        // and the current trackFraction
        // i.e.
        // let s = stepFraction, t = trackFraction
        // Motion of x in time:
        // |-----------------|---------|---------|
        // x00               x0        xi        x
        //
        // where xi is the correct value of x at the required
        // tracking instant.
        //
        // x0 = x00 + s*(x - x00) = s*x + (1 - s)*x00
        //
        // i.e. the motion covered by previous tracking portions
        // within this timestep, and
        //
        // xi = x0 + t*(x - x0)
        //    = t*x + (1 - t)*x0
        //    = t*x + (1 - t)*(s*x + (1 - s)*x00)
        //    = (s + t - s*t)*x + (1 - (s + t - s*t))*x00
        //
        // let m = (s + t - s*t)
        //
        // xi = m*x + (1 - m)*x00 = x00 + m*(x - x00);
        //
        // In the same form as before.

        scalar m =
            p.stepFraction()
          + trackFraction
          - (p.stepFraction()*trackFraction);

        // When the mesh is moving, the velocity field on wall patches
        // will contain the velocity associated with the motion of the
        // mesh, in which case it is interpolated in time using m.
        // For other patches the face velocity will need to be
        // reconstructed from the face centre motion.

        const vector& Cf = mesh.faceCentres()[p.face()];

        vector Cf00 = mesh.faces()[p.face()].centre(mesh.oldPoints());

        if (isA<wallPolyPatch>(pp))
        {
            Up = U00 + m*(U - U00);
        }
        else
        {
            Up = (Cf - Cf00)/this->owner().time().deltaTValue();
        }

        if (mag(dn) > SMALL)
        {
            // Rotational motion, nw requires interpolation and a
            // rotational velocity around face centre correction to Up
            // is required.

            nw = n00 + m*dn;

            // Cf at tracking instant
            vector Cfi = Cf00 + m*(Cf - Cf00);

            // Normal vector cross product
            vector omega = (n00 ^ n);

            scalar magOmega = mag(omega);

            // magOmega = sin(angle between unit normals)
            // Normalise omega vector by magOmega, then multiply by
            // angle/dt to give the correct angular velocity vector.
            omega *=
                Foam::asin(magOmega)
               /(magOmega*this->owner().time().deltaTValue());

            // Project position onto face and calculate this position
            // relative to the face centre.
            vector facePos =
                p.position()
              - ((p.position() - Cfi) & nw)*nw
              - Cfi;

            Up += (omega ^ facePos);
        }

        // No further action is required if the motion is
        // translational only, nw and Up have already been set.
    }
}


template<class CloudType>
void Foam::PatchInteractionModel<CloudType>::info(Ostream& os)
{
    // do nothing
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PatchInteractionModelNew.C"

// ************************************************************************* //
