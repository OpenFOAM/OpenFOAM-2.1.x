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

#include "specieThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo>
inline Foam::specieThermo<Thermo>::specieThermo
(
    const Thermo& sp
)
:
    Thermo(sp)
{}


template<class Thermo>
inline Foam::scalar Foam::specieThermo<Thermo>::T
(
    scalar f,
    scalar T0,
    scalar (specieThermo<Thermo>::*F)(const scalar) const,
    scalar (specieThermo<Thermo>::*dFdT)(const scalar) const,
    scalar (specieThermo<Thermo>::*limit)(const scalar) const
) const
{
    scalar Test = T0;
    scalar Tnew = T0;
    scalar Ttol = T0*tol_;
    int    iter = 0;

    do
    {
        Test = Tnew;
        Tnew =
            (this->*limit)(Test - ((this->*F)(Test) - f)/(this->*dFdT)(Test));

        if (iter++ > maxIter_)
        {
            FatalErrorIn
            (
                "specieThermo<Thermo>::T(scalar f, scalar T0, "
                "scalar (specieThermo<Thermo>::*F)(const scalar) const, "
                "scalar (specieThermo<Thermo>::*dFdT)(const scalar) const"
                ") const"
            )   << "Maximum number of iterations exceeded"
                << abort(FatalError);
        }

    } while (mag(Tnew - Test) > Ttol);

    return Tnew;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::specieThermo<Thermo>::specieThermo
(
    const word& name,
    const specieThermo& st
)
:
    Thermo(name, st)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::specieThermo<Thermo>::cv(const scalar T) const
{
    return this->cp(T) - this->RR;
}


template<class Thermo>
inline Foam::scalar Foam::specieThermo<Thermo>::gamma(const scalar T) const
{
    scalar CP = this->cp(T);
    return CP/(CP - this->RR);
}


template<class Thermo>
inline Foam::scalar Foam::specieThermo<Thermo>::e(const scalar T) const
{
    return this->h(T) - this->RR*(T - this->Tstd);
}


template<class Thermo>
inline Foam::scalar Foam::specieThermo<Thermo>::es(const scalar T) const
{
    return this->hs(T) - this->RR*(T - this->Tstd);
}


template<class Thermo>
inline Foam::scalar Foam::specieThermo<Thermo>::g(const scalar T) const
{
    return this->h(T) - T*this->s(T);
}


template<class Thermo>
inline Foam::scalar Foam::specieThermo<Thermo>::a(const scalar T) const
{
    return this->e(T) - T*this->s(T);
}


template<class Thermo>
inline Foam::scalar Foam::specieThermo<Thermo>::Cp(const scalar T) const
{
    return this->cp(T)/this->W();
}


template<class Thermo>
inline Foam::scalar Foam::specieThermo<Thermo>::Cv(const scalar T) const
{
    return this->cv(T)/this->W();
}


template<class Thermo>
inline Foam::scalar Foam::specieThermo<Thermo>::H(const scalar T) const
{
    return this->h(T)/this->W();
}


template<class Thermo>
inline Foam::scalar Foam::specieThermo<Thermo>::Hs(const scalar T) const
{
    return this->hs(T)/this->W();
}


template<class Thermo>
inline Foam::scalar Foam::specieThermo<Thermo>::Hc() const
{
    return this->hc()/this->W();
}


template<class Thermo>
inline Foam::scalar Foam::specieThermo<Thermo>::S(const scalar T) const
{
    return this->s(T)/this->W();
}


template<class Thermo>
inline Foam::scalar Foam::specieThermo<Thermo>::E(const scalar T) const
{
    return this->e(T)/this->W();
}


template<class Thermo>
inline Foam::scalar Foam::specieThermo<Thermo>::G(const scalar T) const
{
    return this->g(T)/this->W();
}


template<class Thermo>
inline Foam::scalar Foam::specieThermo<Thermo>::A(const scalar T) const
{
    return this->a(T)/this->W();
}


template<class Thermo>
inline Foam::scalar Foam::specieThermo<Thermo>::K(const scalar T) const
{
    if (equal(this->nMoles(), SMALL))
    {
        return 1.0;
    }
    else
    {
        scalar arg = -this->nMoles()*this->g(T)/(this->RR*T);

        if (arg < 600.0)
        {
            return exp(arg);
        }
        else
        {
            return VGREAT;
        }
    }
}


template<class Thermo>
inline Foam::scalar Foam::specieThermo<Thermo>::Kp(const scalar T) const
{
    return K(T);
}


template<class Thermo>
inline Foam::scalar Foam::specieThermo<Thermo>::Kc(const scalar T) const
{
    if (equal(this->nMoles(), SMALL))
    {
        return Kp(T);
    }
    else
    {
        return Kp(T)*pow(this->Pstd/(this->RR*T), this->nMoles());
    }
}


template<class Thermo>
inline Foam::scalar Foam::specieThermo<Thermo>::Kx
(
    const scalar T,
    const scalar p
) const
{
    if (equal(this->nMoles(), SMALL))
    {
        return Kp(T);
    }
    else
    {
        return Kp(T)*pow(this->Pstd/p, this->nMoles());
    }
}


template<class Thermo>
inline Foam::scalar Foam::specieThermo<Thermo>::Kn
(
    const scalar T,
    const scalar p,
    const scalar n
) const
{
    if (equal(this->nMoles(), SMALL))
    {
        return Kp(T);
    }
    else
    {
        return Kp(T)*pow(n*this->Pstd/p, this->nMoles());
    }
}


template<class Thermo>
inline Foam::scalar Foam::specieThermo<Thermo>::TH
(
    const scalar h,
    const scalar T0
) const
{
    return T
    (
        h,
        T0,
        &specieThermo<Thermo>::H,
        &specieThermo<Thermo>::Cp,
        &specieThermo<Thermo>::limit
    );
}


template<class Thermo>
inline Foam::scalar Foam::specieThermo<Thermo>::THs
(
    const scalar hs,
    const scalar T0
) const
{
    return T
    (
        hs,
        T0,
        &specieThermo<Thermo>::Hs,
        &specieThermo<Thermo>::Cp,
        &specieThermo<Thermo>::limit
    );
}


template<class Thermo>
inline Foam::scalar Foam::specieThermo<Thermo>::TE
(
    const scalar e,
    const scalar T0
) const
{
    return T
    (
        e,
        T0,
        &specieThermo<Thermo>::E,
        &specieThermo<Thermo>::Cv,
        &specieThermo<Thermo>::limit
    );
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline void Foam::specieThermo<Thermo>::operator+=
(
    const specieThermo<Thermo>& st
)
{
    Thermo::operator+=(st);
}


template<class Thermo>
inline void Foam::specieThermo<Thermo>::operator-=
(
    const specieThermo<Thermo>& st
)
{
    Thermo::operator-=(st);
}


template<class Thermo>
inline void Foam::specieThermo<Thermo>::operator*=(const scalar s)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::specieThermo<Thermo> Foam::operator+
(
    const specieThermo<Thermo>& st1,
    const specieThermo<Thermo>& st2
)
{
    return specieThermo<Thermo>
    (
        static_cast<const Thermo&>(st1) + static_cast<const Thermo&>(st2)
    );
}


template<class Thermo>
inline Foam::specieThermo<Thermo> Foam::operator-
(
    const specieThermo<Thermo>& st1,
    const specieThermo<Thermo>& st2
)
{
    return specieThermo<Thermo>
    (
        static_cast<const Thermo&>(st1) - static_cast<const Thermo&>(st2)
    );
}


template<class Thermo>
inline Foam::specieThermo<Thermo> Foam::operator*
(
    const scalar s,
    const specieThermo<Thermo>& st
)
{
    return specieThermo<Thermo>
    (
        s*static_cast<const Thermo&>(st)
    );
}


template<class Thermo>
inline Foam::specieThermo<Thermo> Foam::operator==
(
    const specieThermo<Thermo>& st1,
    const specieThermo<Thermo>& st2
)
{
    return st2 - st1;
}


// ************************************************************************* //
