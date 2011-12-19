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

#include "PressureGradientForce.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PressureGradientForce<CloudType>::PressureGradientForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ParticleForce<CloudType>(owner, mesh, dict, typeName, true),
    UName_(this->coeffs().lookup("U")),
    gradUPtr_(NULL)
{}


template<class CloudType>
Foam::PressureGradientForce<CloudType>::PressureGradientForce
(
    const PressureGradientForce& pgf
)
:
    ParticleForce<CloudType>(pgf),
    UName_(pgf.UName_),
    gradUPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PressureGradientForce<CloudType>::~PressureGradientForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::PressureGradientForce<CloudType>::cacheFields(const bool store)
{
    if (store)
    {
        const volVectorField& U = this->mesh().template
            lookupObject<volVectorField>(UName_);
        gradUPtr_ = fvc::grad(U).ptr();
    }
    else
    {
        if (gradUPtr_)
        {
            delete gradUPtr_;
            gradUPtr_ = NULL;
        }
    }
}


template<class CloudType>
Foam::forceSuSp Foam::PressureGradientForce<CloudType>::calcCoupled
(
    const typename CloudType::parcelType& p,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    forceSuSp value(vector::zero, 0.0);

    const volTensorField& gradU = *gradUPtr_;
    value.Su() = mass*p.rhoc()/p.rho()*(p.U() & gradU[p.cell()]);

    return value;
}


// ************************************************************************* //
