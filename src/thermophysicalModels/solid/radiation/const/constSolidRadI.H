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

#include "constSolidRad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class thermo>
inline constSolidRad<thermo>::constSolidRad
(
    const thermo& t,
    const scalar kappa,
    const scalar sigmaS,
    const scalar emissivity
)
:
    thermo(t),
    kappa_(kappa),
    sigmaS_(sigmaS),
    emissivity_(emissivity)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


template<class thermo>
inline constSolidRad<thermo>::constSolidRad
(
    const word& name,
    const constSolidRad& pg
)
:
    thermo(name, pg),
    kappa_(pg.kappa_),
    sigmaS_(pg.sigmaS_),
    emissivity_(pg.emissivity_)
{}


template<class thermo>
inline autoPtr<constSolidRad<thermo> > constSolidRad<thermo>::clone() const
{
    return autoPtr<constSolidRad<thermo> >(new constSolidRad<thermo>(*this));
}


template<class thermo>
inline autoPtr<constSolidRad<thermo> > constSolidRad<thermo>::New(Istream& is)
{
    return autoPtr<constSolidRad<thermo> >(new constSolidRad<thermo>(is));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class thermo>
inline scalar constSolidRad<thermo>::kappa(scalar) const
{
    return kappa_;
}


template<class thermo>
inline scalar constSolidRad<thermo>::sigmaS(scalar) const
{
    return sigmaS_;
}


template<class thermo>
inline scalar constSolidRad<thermo>::emissivity(scalar) const
{
    return emissivity_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
