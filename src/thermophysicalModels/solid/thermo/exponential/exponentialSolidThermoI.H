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

#include "exponentialSolidThermo.H"
#include "specie.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


template<class rhoType>
inline void Foam::exponentialSolidThermo<rhoType>::checkT(const scalar T) const
{
    if (T < 0.0)
    {
        FatalErrorIn
        (
            "exponentialSolidThermo<rhoType>::checkT(const scalar T) const"
        )   << "attempt to use exponentialSolidThermo<rhoType>"
               " out of temperature range "
            << T
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class rhoType>
inline Foam::exponentialSolidThermo<rhoType>::exponentialSolidThermo
(
    const word& name,
    const exponentialSolidThermo& jt
)
:
    rhoType(name, jt),
    c0_(jt.c0_),
    n0_(jt.n0_),
    Tref_(jt.Tref_),
    Hf_(jt.Hf_)
{}


template<class rhoType>
inline Foam::exponentialSolidThermo<rhoType>::exponentialSolidThermo
(
    const rhoType& st,
    const scalar c0,
    const scalar n0,
    const scalar Tref,
    const scalar Hf
)
:
    rhoType(st),
    c0_(c0),
    n0_(n0),
    Tref_(Tref),
    Hf_(Hf)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class rhoType>
inline Foam::scalar Foam::exponentialSolidThermo<rhoType>::integrateCp
(
    const scalar T
) const
{
    return
    (
        c0_*pow(T, n0_ + 1.0)
       /(pow(Tref_, n0_)*(n0_ + 1.0))
    );
}


template<class rhoType>
inline Foam::scalar Foam::exponentialSolidThermo<rhoType>::Cp
(
    const scalar T
) const
{
    return c0_*pow(T/Tref_, n0_);
}


template<class rhoType>
inline Foam::scalar Foam::exponentialSolidThermo<rhoType>::h
(
    const scalar T
) const
{
    scalar hOffset = integrateCp(specie::Tstd);

    return
    (
        integrateCp(T) + Hf_ - hOffset
    );
}


template<class rhoType>
inline Foam::scalar Foam::exponentialSolidThermo<rhoType>::hs
(
    const scalar T
) const
{
    return h(T) - Hf_;
}


template<class rhoType>
inline Foam::scalar Foam::exponentialSolidThermo<rhoType>::hf() const
{
    return Hf_;
}


// ************************************************************************* //
