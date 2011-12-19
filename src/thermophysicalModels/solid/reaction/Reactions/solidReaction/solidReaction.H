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
    Foam::solidReaction

Description

    Read solid reactions of the type S1 = S2 + G1

SourceFiles
    solidReactionI.H
    solidReaction.C

\*---------------------------------------------------------------------------*/

#ifndef solidReaction_H
#define solidReaction_H

#include "speciesTable.H"
#include "scalarField.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class solidReaction;

inline Ostream& operator<<(Ostream&, const solidReaction&);

/*---------------------------------------------------------------------------*\
                        Class solidReaction Declaration
\*---------------------------------------------------------------------------*/

class solidReaction
{

private:

    // Private data

        //- List of solid names present in reaction system
        const speciesTable& components_;

        //- List of gas species present in reaction system
        speciesTable pyrolisisGases_;

        //- Solid components index for the left-hand-side of the reaction
        List<label> slhs_;

        //- Solid components index for the right-hand-side of the reaction
        List<label> srhs_;

        //- Specie index for the right-hand-side of the reaction
        List<label> grhs_;


    // Private Member Functions

        //- Set rhs and lhs of the reaction
        void setLRhs(Istream&);

        //- Look for the component index in the reaction
        label componentIndex(bool& isGas, Istream& is);

        //- Disallow default bitwise assignment
        void operator=(const solidReaction&);


public:

    //- Runtime type information
    TypeName("Reaction");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            solidReaction,
            Istream,
            (
                const speciesTable& components,
                Istream& is,
                const speciesTable& pyrolysisGases
            ),
            (components, is, pyrolysisGases)
        );


    // Public classes

        //- Class used for the read-construction of PtrLists of reaction
        class iNew
        {
            const speciesTable& components_;
            speciesTable pyrolisisGases_;

        public:

            iNew
            (
                const speciesTable& components,
                Istream& pyrolisisGases
            )
            :
                components_(components),
                pyrolisisGases_(pyrolisisGases)
            {}

            autoPtr<solidReaction> operator()(Istream& is) const
            {
                return autoPtr<solidReaction>
                (
                    solidReaction::New
                    (
                        components_,
                        is,
                        pyrolisisGases_
                    )
                );
            }
        };


    // Constructors

        //- Construct from components
        solidReaction
        (
            const speciesTable& components,
            const speciesTable& pyrolisisGases,
            const List<label>& slhs,
            const List<label>& srhs,
            const List<label>& grhs
        );

        //- Construct as copy given new speciesTable
        solidReaction
        (
            const solidReaction&,
            const speciesTable& components,
            const speciesTable& pyrolisisGases
        );

        //- Construct from Istream
        solidReaction
        (
            const speciesTable& components,
            Istream& is,
            const speciesTable& pyrolisisGases
        );

        //- Construct and return a clone
        virtual autoPtr<solidReaction > clone() const
        {
            return autoPtr<solidReaction >
            (
                new solidReaction(*this)
            );
        }


    // Selectors

        //- Return a pointer to a new patchField created on freestore from input
        static autoPtr<solidReaction > New
        (
            const speciesTable& components,
            Istream&,
            const speciesTable& pyrolisisGases
        );


    //- Destructor
    virtual ~solidReaction()
    {}


    // Member Functions

        // Access

            inline const List<label>& slhs() const;
            inline const List<label>& srhs() const;
            inline const List<label>& grhs() const;

            inline const speciesTable& pyrolisisGases() const;


        // solidReaction rate coefficients

            virtual scalar kf
            (
                const scalar T,
                const scalar p,
                const scalarField& c
            ) const;

            virtual scalar nReact() const;


        //- Write
        virtual void write(Ostream&) const;


    // Ostream Operator

        friend Ostream& operator<<
        (
            Ostream&,
            const solidReaction&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "solidReactionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
