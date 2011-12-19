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

#include "IrreversibleSolidReaction.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionRate>
Foam::IrreversibleSolidReaction<ReactionRate>::IrreversibleSolidReaction
(
    const solidReaction& reaction,
    const ReactionRate& k,
    const scalar nReact
)
:
    solidReaction(reaction),
    k_(k),
    nReact_(nReact)
{}


template<class ReactionRate>
Foam::IrreversibleSolidReaction<ReactionRate>::IrreversibleSolidReaction
(
    const speciesTable& components,
    Istream& is,
    const speciesTable& pyrolysisGases
)
:
    solidReaction(components, is, pyrolysisGases),
    k_(components, is),
    nReact_(readScalar(is))
{
    is.readEnd("solidArrheniusReactionRate(Istream&)");
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ReactionRate>
Foam::scalar Foam::IrreversibleSolidReaction<ReactionRate>::kf
(
    const scalar T,
    const scalar p,
    const scalarField& c
) const
{
    return k_(T, p, c);
}


template<class ReactionRate>
Foam::scalar Foam::IrreversibleSolidReaction<ReactionRate>::nReact() const
{
    return nReact_;
}


template<class ReactionRate>
void Foam::IrreversibleSolidReaction<ReactionRate>::write
(
    Ostream& os
) const
{
    solidReaction::write(os);
    os  << token::SPACE << "Reaction order: " << nReact_ << nl << k_;
}


// ************************************************************************* //
