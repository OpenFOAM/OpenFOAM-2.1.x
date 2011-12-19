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

Class
    Foam::IrreversibleSolidReaction

Description
    Simple extension of Reaction to handle reversible reactions

SourceFiles
    IrreversibleSolidReaction.C

\*---------------------------------------------------------------------------*/

#ifndef IrreversibleSolidReaction_H
#define IrreversibleSolidReaction_H

#include "solidReaction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class IrreversibleSolidReaction Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionRate>
class IrreversibleSolidReaction
:
    public solidReaction
{
    // Private data

        // Reaction rate
        ReactionRate k_;

        // Reaction order
        scalar nReact_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=
        (
            const IrreversibleSolidReaction<ReactionRate>&
        );


public:

    //- Runtime type information
    TypeName("irreversible");


    // Constructors

        //- Construct from components
        IrreversibleSolidReaction
        (
            const solidReaction& reaction,
            const ReactionRate& reactionRate,
            const scalar nReact
        );


        //- Construct from Istream
        IrreversibleSolidReaction
        (
            const speciesTable& components,
            Istream& is,
            const speciesTable& pyrolysisGases
        );


    //- Destructor
    virtual ~IrreversibleSolidReaction()
    {}


    // Member Functions

        // IrreversibleSolidReaction rate coefficients

            //- Forward rate constant
            virtual scalar kf
            (
                const scalar T,
                const scalar p,
                const scalarField& c
            ) const;


            //- Reaction order
            virtual scalar nReact() const;


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "IrreversibleSolidReaction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
