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

#include "NonEquilibriumReversibleReaction.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo, class ReactionRate>
Foam::NonEquilibriumReversibleReaction<ReactionThermo, ReactionRate>::
NonEquilibriumReversibleReaction
(
    const Reaction<ReactionThermo>& reaction,
    const ReactionRate& forwardReactionRate,
    const ReactionRate& reverseReactionRate
)
:
    Reaction<ReactionThermo>(reaction),
    fk_(forwardReactionRate),
    rk_(reverseReactionRate)
{}



template<class ReactionThermo, class ReactionRate>
Foam::NonEquilibriumReversibleReaction<ReactionThermo, ReactionRate>::
NonEquilibriumReversibleReaction
(
    const speciesTable& species,
    const HashPtrTable<ReactionThermo>& thermoDatabase,
    Istream& is
)
:
    Reaction<ReactionThermo>(species, thermoDatabase, is),
    fk_(species, is),
    rk_(species, is)
{}


template<class ReactionThermo, class ReactionRate>
Foam::NonEquilibriumReversibleReaction<ReactionThermo, ReactionRate>::
NonEquilibriumReversibleReaction
(
    const speciesTable& species,
    const HashPtrTable<ReactionThermo>& thermoDatabase,
    const dictionary& dict
)
:
    Reaction<ReactionThermo>(species, thermoDatabase, dict),
    fk_(species, dict.subDict("forward")),
    rk_(species, dict.subDict("reverse"))
{}


template<class ReactionThermo, class ReactionRate>
Foam::NonEquilibriumReversibleReaction<ReactionThermo, ReactionRate>::
NonEquilibriumReversibleReaction
(
    const NonEquilibriumReversibleReaction<ReactionThermo, ReactionRate>& nerr,
    const speciesTable& species
)
:
    Reaction<ReactionThermo>(nerr, species),
    fk_(nerr.fk_),
    rk_(nerr.rk_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ReactionThermo, class ReactionRate>
Foam::scalar
Foam::NonEquilibriumReversibleReaction<ReactionThermo, ReactionRate>::kf
(
    const scalar T,
    const scalar p,
    const scalarField& c
) const
{
    return fk_(T, p, c);
}


template<class ReactionThermo, class ReactionRate>
Foam::scalar
Foam::NonEquilibriumReversibleReaction<ReactionThermo, ReactionRate>::kr
(
    const scalar,
    const scalar T,
    const scalar p,
    const scalarField& c
) const
{
    return rk_(T, p, c);
}


template<class ReactionThermo, class ReactionRate>
Foam::scalar
Foam::NonEquilibriumReversibleReaction<ReactionThermo, ReactionRate>::kr
(
    const scalar T,
    const scalar p,
    const scalarField& c
) const
{
    return rk_(T, p, c);
}


template<class ReactionThermo, class ReactionRate>
void Foam::NonEquilibriumReversibleReaction<ReactionThermo, ReactionRate>::write
(
    Ostream& os
) const
{
    Reaction<ReactionThermo>::write(os);

    os  << indent << "forward" << nl;
    os  << indent << token::BEGIN_BLOCK << nl;
    os  << incrIndent;
    fk_.write(os);
    os  << decrIndent;
    os  << indent << token::END_BLOCK << nl;

    os  << indent << "reverse" << nl;
    os  << indent << token::BEGIN_BLOCK << nl;
    os  << incrIndent;
    rk_.write(os);
    os  << decrIndent;
    os  << indent << token::END_BLOCK << nl;
}


// ************************************************************************* //
