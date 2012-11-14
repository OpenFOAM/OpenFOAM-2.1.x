/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "TableBase.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::TableBase<Type>::TableBase(const word& name, const dictionary& dict)
:
    name_(name),
    boundsHandling_
    (
        wordToBoundsHandling
        (
            dict.lookupOrDefault<word>("outOfBounds", "clamp")
        )
    ),
    table_()
{}


template<class Type>
Foam::TableBase<Type>::TableBase(const TableBase<Type>& tbl)
:
    name_(tbl.name_),
    boundsHandling_(tbl.boundsHandling_),
    table_(tbl.table_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::TableBase<Type>::~TableBase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::word Foam::TableBase<Type>::boundsHandlingToWord
(
     const boundsHandling& bound
) const
{
    word enumName("warn");

    switch (bound)
    {
        case ERROR:
        {
            enumName = "error";
            break;
        }
        case WARN:
        {
            enumName = "warn";
            break;
        }
        case CLAMP:
        {
            enumName = "clamp";
            break;
        }
        case REPEAT:
        {
            enumName = "repeat";
            break;
        }
    }

    return enumName;
}


template<class Type>
typename Foam::TableBase<Type>::boundsHandling
Foam::TableBase<Type>::wordToBoundsHandling
(
    const word& bound
) const
{
    if (bound == "error")
    {
        return ERROR;
    }
    else if (bound == "warn")
    {
        return WARN;
    }
    else if (bound == "clamp")
    {
        return CLAMP;
    }
    else if (bound == "repeat")
    {
        return REPEAT;
    }
    else
    {
        WarningIn("Foam::TableBase<Type>::wordToBoundsHandling(const word&)")
            << "bad outOfBounds specifier " << bound << " using 'warn'"
            << endl;

        return WARN;
    }
}


template<class Type>
typename Foam::TableBase<Type>::boundsHandling
Foam::TableBase<Type>::outOfBounds
(
    const boundsHandling& bound
)
{
    boundsHandling prev = boundsHandling_;
    boundsHandling_ = bound;
    return prev;
}


template<class Type>
void Foam::TableBase<Type>::check() const
{
    if (!table_.size())
    {
        FatalErrorIn("Foam::TableBase<Type>::check() const")
            << "Table for entry " << this->name_ << " is invalid (empty)"
            << nl << exit(FatalError);
    }

    label n = table_.size();
    scalar prevValue = table_[0].first();

    for (label i = 1; i < n; ++i)
    {
        const scalar currValue = table_[i].first();

        // avoid duplicate values (divide-by-zero error)
        if (currValue <= prevValue)
        {
            FatalErrorIn("Foam::TableBase<Type>::check() const")
                << "out-of-order value: " << currValue << " at index " << i
                << exit(FatalError);
        }
        prevValue = currValue;
    }
}


template<class Type>
bool Foam::TableBase<Type>::checkMinBounds
(
    const scalar x,
    scalar& xDash
) const
{
    if (x < table_[0].first())
    {
        switch (boundsHandling_)
        {
            case ERROR:
            {
                FatalErrorIn
                (
                    "bool Foam::TableBase<Type>::checkMinBounds"
                    "("
                        "const scalar, "
                        "scalar&"
                    ") const"
                )   << "value (" << x << ") underflow"
                    << exit(FatalError);
                break;
            }
            case WARN:
            {
                WarningIn
                (
                    "bool Foam::TableBase<Type>::checkMinBounds"
                    "("
                        "const scalar, "
                        "scalar&"
                    ") const"
                )   << "value (" << x << ") underflow" << nl
                    << endl;

                // fall-through to 'CLAMP'
            }
            case CLAMP:
            {
                xDash = table_[0].first();
                return true;
                break;
            }
            case REPEAT:
            {
                // adjust x to >= minX
                scalar span = table_.last().first() - table_[0].first();
                xDash = fmod(x - table_[0].first(), span) + table_[0].first();
                break;
            }
        }
    }
    else
    {
        xDash = x;
    }

    return false;
}


template<class Type>
bool Foam::TableBase<Type>::checkMaxBounds
(
    const scalar x,
    scalar& xDash
) const
{
    if (x > table_.last().first())
    {
        switch (boundsHandling_)
        {
            case ERROR:
            {
                FatalErrorIn
                (
                    "bool Foam::TableBase<Type>::checkMaxBounds"
                    "("
                        "const scalar, "
                        "scalar&"
                    ") const"
                )   << "value (" << x << ") overflow"
                    << exit(FatalError);
                break;
            }
            case WARN:
            {
                WarningIn
                (
                    "bool Foam::TableBase<Type>::checkMaxBounds"
                    "("
                        "const scalar, "
                        "scalar&"
                    ") const"
                )   << "value (" << x << ") overflow" << nl
                    << endl;

                // fall-through to 'CLAMP'
            }
            case CLAMP:
            {
                xDash = table_.last().first();
                return true;
                break;
            }
            case REPEAT:
            {
                // adjust x to >= minX
                scalar span = table_.last().first() - table_[0].first();
                xDash = fmod(x - table_[0].first(), span) + table_[0].first();
                break;
            }
        }
    }
    else
    {
        xDash = x;
    }

    return false;
}


template<class Type>
void Foam::TableBase<Type>::convertTimeBase(const Time& t)
{
    forAll(table_, i)
    {
        scalar value = table_[i].first();
        table_[i].first() = t.userTimeToTime(value);
    }
}


template<class Type>
Type Foam::TableBase<Type>::value(const scalar x) const
{
    scalar xDash = x;

    if (checkMinBounds(x, xDash))
    {
        return table_[0].second();
    }

    if (checkMaxBounds(xDash, xDash))
    {
        return table_.last().second();
    }

    // Find i such that x(i) < xDash < x(i+1)
    label i = 0;
    while ((i+1 < table_.size()) && (table_[i+1].first() < xDash))
    {
        i++;
    }

    if (i+1 == table_.size())
    {
        // Reached end of table. This can happen for tables of size 1.
        return table_[i].second();
    }

    // Linear interpolation to find value
    return Type
    (
        (xDash - table_[i].first())/(table_[i+1].first() - table_[i].first())
      * (table_[i+1].second() - table_[i].second())
      + table_[i].second()
    );
}


template<class Type>
Type Foam::TableBase<Type>::integrate(const scalar x1, const scalar x2) const
{
    // Initialise return value
    Type sum = pTraits<Type>::zero;

    // Return zero if out of bounds
    if ((x1 > table_.last().first()) || (x2 < table_[0].first()))
    {
        return sum;
    }

    // Find next index greater than x1
    label id1 = 0;
    while ((table_[id1].first() < x1) && (id1 < table_.size()))
    {
        id1++;
    }

    // Find next index less than x2
    label id2 = table_.size() - 1;
    while ((table_[id2].first() > x2) && (id2 >= 1))
    {
        id2--;
    }

    if ((id1 - id2) == 1)
    {
        // x1 and x2 lie within 1 interval
        sum = 0.5*(value(x1) + value(x2))*(x2 - x1);
    }
    else
    {
        // x1 and x2 cross multiple intervals

        // Integrate table body
        for (label i=id1; i<id2; i++)
        {
            sum +=
                (table_[i].second() + table_[i+1].second())
              * (table_[i+1].first() - table_[i].first());
        }
        sum *= 0.5;

        // Add table ends (partial segments)
        if (id1 > 0)
        {
            sum += 0.5
              * (value(x1) + table_[id1].second())
              * (table_[id1].first() - x1);
        }
        if (id2 < table_.size() - 1)
        {
            sum += 0.5
              * (table_[id2].second() + value(x2))
              * (x2 - table_[id2].first());
        }
    }

    return sum;
}


// * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * * * //

#include "TableBaseIO.C"


// ************************************************************************* //
