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
    Foam::reactingSolidMixture

Description
    Foam::reactingSolidMixture

SourceFiles
    reactingSolidMixture.C

\*---------------------------------------------------------------------------*/

#ifndef reactingSolidMixture_H
#define reactingSolidMixture_H

#include "multiComponentSolidMixture.H"
#include "solidReaction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class reactingSolidMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoSolidType>
class reactingSolidMixture
:
    public multiComponentSolidMixture<ThermoSolidType>,
    public PtrList<solidReaction>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        reactingSolidMixture(const reactingSolidMixture&);

        //- Disallow default bitwise assignment
        void operator=(const reactingSolidMixture&);


public:

    //- The type of thermo package this mixture is instantiated for
    typedef ThermoSolidType thermoType;


    // Constructors

        //- Construct from dictionary and mesh
        reactingSolidMixture(const dictionary&, const fvMesh&);


    //- Destructor
    virtual ~reactingSolidMixture()
    {}


    // Member functions

        //- Read dictionary
        void read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "reactingSolidMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
