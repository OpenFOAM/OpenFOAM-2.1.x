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

Class
    Foam::isotropicKSolidThermo

Description
    Directional conductivity + table interpolation.

SourceFiles
    isotropicKSolidThermo.C

\*---------------------------------------------------------------------------*/

#ifndef directionalSolidThermo_H
#define directionalSolidThermo_H

#include "interpolatedSolidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class isotropicKSolidThermo Declaration
\*---------------------------------------------------------------------------*/

class isotropicKSolidThermo
:
    public interpolatedSolidThermo
{
    // Private data

        //- Thermal conductivity [W/(m.K)]
        volScalarField K_;

        //- Thermal conductivity vector
        Field<scalar> KValues_;


public:

    //- Runtime type information
    TypeName("isotropicKSolidThermo");


    // Constructors

        //- Construct from mesh
        isotropicKSolidThermo(const fvMesh& mesh);

        //- Construct from mesh and dicionary
        isotropicKSolidThermo(const fvMesh& mesh, const dictionary& dict);



    //- Destructor
    virtual ~isotropicKSolidThermo();


    // Member functions

        //- Update properties
        virtual void correct();

        // Access functions

            //- Constant access to K
            virtual tmp<volScalarField> K() const;

            //- Thermal conductivity [W/m/K]
            virtual tmp<volSymmTensorField> directionalK() const
            {
                notImplemented("directionalKSolidThermo::directionalK() const");
                return tmp<volSymmTensorField>(NULL);
            }

        // Per patch calculation

            //- Thermal conductivity [W//m/K]
            virtual tmp<scalarField> K(const label patchI) const;

            //- Thermal conductivity [W//m/K]
            virtual tmp<symmTensorField> directionalK(const label) const
            {
                notImplemented
                (
                    "directionalKSolidThermo::directionalK(const label) const"
                );
                return tmp<symmTensorField>(NULL);
            }


    // I-O

        //- Read solidThermophysicalProperties dictionary
        virtual bool read();

        //- Write properties
        virtual bool writeData(Ostream& os) const;

        //- Ostream Operator
        friend Ostream& operator<<
        (
            Ostream& os,
            const isotropicKSolidThermo& s
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
