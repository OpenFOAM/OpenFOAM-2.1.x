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

#include "icoPolynomial.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<int PolySize>
inline Foam::icoPolynomial<PolySize>::icoPolynomial
(
    const specie& sp,
    const Polynomial<PolySize>& rhoCoeffs
)
:
    specie(sp),
    rhoCoeffs_(rhoCoeffs)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<int PolySize>
inline Foam::icoPolynomial<PolySize>::icoPolynomial
(
    const icoPolynomial<PolySize>& ip
)
:
    specie(ip),
    rhoCoeffs_(ip.rhoCoeffs_)
{}


template<int PolySize>
inline Foam::icoPolynomial<PolySize>::icoPolynomial
(
    const word& name,
    const icoPolynomial<PolySize>& ip
)
:
    specie(name, ip),
    rhoCoeffs_(ip.rhoCoeffs_)
{}


template<int PolySize>
inline Foam::autoPtr<Foam::icoPolynomial<PolySize> >
Foam::icoPolynomial<PolySize>::clone() const
{
    return autoPtr<icoPolynomial<PolySize> >
    (
        new icoPolynomial<PolySize>(*this)
    );
}


template<int PolySize>
inline Foam::autoPtr<Foam::icoPolynomial<PolySize> >
Foam::icoPolynomial<PolySize>::New(Istream& is)
{
    return autoPtr<icoPolynomial<PolySize> >(new icoPolynomial<PolySize>(is));
}


template<int PolySize>
inline Foam::autoPtr<Foam::icoPolynomial<PolySize> >
Foam::icoPolynomial<PolySize>::New(const dictionary& dict)
{
    return autoPtr<icoPolynomial<PolySize> >
    (
        new icoPolynomial<PolySize>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<int PolySize>
inline Foam::scalar Foam::icoPolynomial<PolySize>::rho(scalar, scalar T) const
{
    return rhoCoeffs_.value(T)/this->W();
}


template<int PolySize>
inline Foam::scalar Foam::icoPolynomial<PolySize>::psi(scalar, scalar) const
{
    return 0.0;
}


template<int PolySize>
inline Foam::scalar Foam::icoPolynomial<PolySize>::Z(scalar, scalar) const
{
    return 0.0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<int PolySize>
inline Foam::icoPolynomial<PolySize>& Foam::icoPolynomial<PolySize>::operator=
(
    const icoPolynomial<PolySize>& ip
)
{
    specie::operator=(ip);

    rhoCoeffs_ = ip.rhoCoeffs_;

    return *this;
}


template<int PolySize>
inline void Foam::icoPolynomial<PolySize>::operator+=
(
    const icoPolynomial<PolySize>& ip
)
{
    scalar molr1 = this->nMoles();

    specie::operator+=(ip);

    molr1 /= this->nMoles();
    scalar molr2 = ip.nMoles()/this->nMoles();

    rhoCoeffs_ = molr1*rhoCoeffs_ + molr2*ip.rhoCoeffs_;
}


template<int PolySize>
inline void Foam::icoPolynomial<PolySize>::operator-=
(
    const icoPolynomial<PolySize>& ip
)
{
    scalar molr1 = this->nMoles();

    specie::operator-=(ip);

    molr1 /= this->nMoles();
    scalar molr2 = ip.nMoles()/this->nMoles();

    rhoCoeffs_ = molr1*rhoCoeffs_ - molr2*ip.rhoCoeffs_;
}


template<int PolySize>
inline void Foam::icoPolynomial<PolySize>::operator*=(const scalar s)
{
    specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<int PolySize>
Foam::icoPolynomial<PolySize> Foam::operator+
(
    const icoPolynomial<PolySize>& ip1,
    const icoPolynomial<PolySize>& ip2
)
{
    scalar nMoles = ip1.nMoles() + ip2.nMoles();
    scalar molr1 = ip1.nMoles()/nMoles;
    scalar molr2 = ip2.nMoles()/nMoles;

    return icoPolynomial<PolySize>
    (
        static_cast<const specie&>(ip1)
      + static_cast<const specie&>(ip2),
        molr1*ip1.rhoCoeffs_ + molr2*ip2.rhoCoeffs_
    );
}


template<int PolySize>
Foam::icoPolynomial<PolySize> Foam::operator-
(
    const icoPolynomial<PolySize>& ip1,
    const icoPolynomial<PolySize>& ip2
)
{
    scalar nMoles = ip1.nMoles() + ip2.nMoles();
    scalar molr1 = ip1.nMoles()/nMoles;
    scalar molr2 = ip2.nMoles()/nMoles;

    return icoPolynomial<PolySize>
    (
        static_cast<const specie&>(ip1)
      - static_cast<const specie&>(ip2),
        molr1*ip1.rhoCoeffs_ - molr2*ip2.rhoCoeffs_
    );
}


template<int PolySize>
Foam::icoPolynomial<PolySize> Foam::operator*
(
    const scalar s,
    const icoPolynomial<PolySize>& ip
)
{
    return icoPolynomial<PolySize>
    (
        s*static_cast<const specie&>(ip),
        ip.rhoCoeffs_
    );
}


template<int PolySize>
Foam::icoPolynomial<PolySize> Foam::operator==
(
    const icoPolynomial<PolySize>& ip1,
    const icoPolynomial<PolySize>& ip2
)
{
    return ip2 - ip1;
}


// ************************************************************************* //
