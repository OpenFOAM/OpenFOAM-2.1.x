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

#include "Polynomial.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<int PolySize>
Foam::Polynomial<PolySize>::Polynomial()
:
    VectorSpace<Polynomial<PolySize>, scalar, PolySize>(),
    logActive_(false),
    logCoeff_(0.0)
{
    for (int i = 0; i < PolySize; ++i)
    {
        this->v_[i] = 0.0;
    }
}


template<int PolySize>
Foam::Polynomial<PolySize>::Polynomial
(
    const Polynomial<PolySize>& poly
)
:
    VectorSpace<Polynomial<PolySize>, scalar, PolySize>(poly),
    logActive_(poly.logActive_),
    logCoeff_(poly.logCoeff_)
{}


template<int PolySize>
Foam::Polynomial<PolySize>::Polynomial(const scalar coeffs[PolySize])
:
    VectorSpace<Polynomial<PolySize>, scalar, PolySize>(),
    logActive_(false),
    logCoeff_(0.0)
{
    for (int i=0; i<PolySize; i++)
    {
        this->v_[i] = coeffs[i];
    }
}


template<int PolySize>
Foam::Polynomial<PolySize>::Polynomial(const UList<scalar>& coeffs)
:
    VectorSpace<Polynomial<PolySize>, scalar, PolySize>(),
    logActive_(false),
    logCoeff_(0.0)
{
    if (coeffs.size() != PolySize)
    {
        FatalErrorIn
        (
            "Polynomial<PolySize>::Polynomial(const UList<scalar>&)"
        )   << "Size mismatch: Needed " << PolySize
            << " but given " << coeffs.size()
            << nl << exit(FatalError);
    }

    for (int i = 0; i < PolySize; ++i)
    {
        this->v_[i] = coeffs[i];
    }
}


// template<int PolySize>
// Foam::Polynomial<PolySize>::Polynomial(const polynomialFunction& poly)
// :
//     VectorSpace<Polynomial<PolySize>, scalar, PolySize>(),
//     logActive_(poly.logActive()),
//     logCoeff_(poly.logCoeff())
// {
//     if (poly.size() != PolySize)
//     {
//         FatalErrorIn
//         (
//             "Polynomial<PolySize>::Polynomial(const polynomialFunction&)"
//         )   << "Size mismatch: Needed " << PolySize
//             << " but given " << poly.size()
//             << nl << exit(FatalError);
//     }
//
//     for (int i = 0; i < PolySize; ++i)
//     {
//         this->v_[i] = poly[i];
//     }
// }


template<int PolySize>
Foam::Polynomial<PolySize>::Polynomial(Istream& is)
:
    VectorSpace<Polynomial<PolySize>, scalar, PolySize>(is),
    logActive_(false),
    logCoeff_(0.0)
{}


template<int PolySize>
Foam::Polynomial<PolySize>::Polynomial(const word& name, Istream& is)
:
    VectorSpace<Polynomial<PolySize>, scalar, PolySize>(),
    logActive_(false),
    logCoeff_(0.0)
{
    word isName(is);

    if (isName != name)
    {
        FatalErrorIn
        (
            "Polynomial<PolySize>::Polynomial(const word&, Istream&)"
        )   << "Expected polynomial name " << name << " but read " << isName
            << nl << exit(FatalError);
    }

    VectorSpace<Polynomial<PolySize>, scalar, PolySize>::
        operator=(VectorSpace<Polynomial<PolySize>, scalar, PolySize>(is));

    if (this->size() == 0)
    {
        FatalErrorIn
        (
            "Polynomial<PolySize>::Polynomial(const word&, Istream&)"
        )   << "Polynomial coefficients for entry " << isName
            << " are invalid (empty)" << nl << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<int PolySize>
bool Foam::Polynomial<PolySize>::logActive() const
{
    return logActive_;
}


template<int PolySize>
Foam::scalar Foam::Polynomial<PolySize>::logCoeff() const
{
    return logCoeff_;
}


template<int PolySize>
Foam::scalar Foam::Polynomial<PolySize>::value(const scalar x) const
{
    scalar val = this->v_[0];

    // avoid costly pow() in calculation
    scalar powX = x;
    for (label i=1; i<PolySize; ++i)
    {
        val += this->v_[i]*powX;
        powX *= x;
    }

    if (logActive_)
    {
        val += logCoeff_*log(x);
    }

    return val;
}


template<int PolySize>
Foam::scalar Foam::Polynomial<PolySize>::integrate
(
    const scalar x1,
    const scalar x2
) const
{
    if (logActive_)
    {
        FatalErrorIn
        (
            "scalar Polynomial<PolySize>::integrate"
            "("
                "const scalar, "
                "const scalar"
            ") const"
        )   << "Cannot integrate polynomial with logarithmic coefficients"
            << nl << abort(FatalError);
    }


    // avoid costly pow() in calculation
    scalar powX1 = x1;
    scalar powX2 = x2;

    scalar val = this->v_[0]*(powX2 - powX1);
    for (label i=1; i<PolySize; ++i)
    {
        val += this->v_[i]/(i + 1) * (powX2 - powX1);
        powX1 *= x1;
        powX2 *= x2;
    }

    return val;
}


template<int PolySize>
typename Foam::Polynomial<PolySize>::intPolyType
Foam::Polynomial<PolySize>::integral(const scalar intConstant) const
{
    intPolyType newCoeffs;

    newCoeffs[0] = intConstant;
    forAll(*this, i)
    {
        newCoeffs[i+1] = this->v_[i]/(i + 1);
    }

    return newCoeffs;
}


template<int PolySize>
typename Foam::Polynomial<PolySize>::polyType
Foam::Polynomial<PolySize>::integralMinus1(const scalar intConstant) const
{
    polyType newCoeffs;

    if (this->v_[0] > VSMALL)
    {
        newCoeffs.logActive_ = true;
        newCoeffs.logCoeff_ = this->v_[0];
    }

    newCoeffs[0] = intConstant;
    for (label i=1; i<PolySize; ++i)
    {
        newCoeffs[i] = this->v_[i]/i;
    }

    return newCoeffs;
}


// ************************************************************************* //
