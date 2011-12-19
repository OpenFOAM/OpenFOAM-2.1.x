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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class rhoType>
inline Foam::constSolidThermo<rhoType>::constSolidThermo
(
    const word& name,
    const constSolidThermo& ct
)
:
    rhoType(name, ct),
    Cp_(ct.Cp_),
    Hf_(ct.Hf_)
{}


template<class rhoType>
inline Foam::constSolidThermo<rhoType>::constSolidThermo
(
    const rhoType& st,
    const scalar Cp,
    const scalar hf
)
:
    rhoType(st),
    Cp_(Cp),
    Hf_(hf)
{}


template<class rhoType>
inline Foam::autoPtr<Foam::constSolidThermo<rhoType> >
Foam::constSolidThermo<rhoType>::clone() const
{
    return autoPtr<constSolidThermo<rhoType> >
    (
        new constSolidThermo<rhoType>(*this)
    );
}


template<class rhoType>
inline Foam::autoPtr<Foam::constSolidThermo<rhoType> >
Foam::constSolidThermo<rhoType>::New(dictionary& is)
{
    return autoPtr<constSolidThermo<rhoType> >
    (
        new constSolidThermo<rhoType>(is)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class rhoType>
inline Foam::scalar Foam::constSolidThermo<rhoType>::Cp
(
    const scalar
) const
{
    return Cp_;
}


template<class rhoType>
inline Foam::scalar Foam::constSolidThermo<rhoType>::h
(
    const scalar T
) const
{
    scalar hOffset = Cp_*(specie::Tstd);

    return (Cp_*T + Hf_ - hOffset);
}


template<class rhoType>
inline Foam::scalar Foam::constSolidThermo<rhoType>::hs
(
    const scalar T
) const
{
    return h(T) - Hf_;
}


template<class rhoType>
inline Foam::scalar Foam::constSolidThermo<rhoType>::hf() const
{
    return Hf_;
}


// ************************************************************************* //
