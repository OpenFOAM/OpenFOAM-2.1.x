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
    Foam::directionalKSolidThermo

Description
    Directional conductivity + table interpolation.

SourceFiles
    directionalKSolidThermo.C

\*---------------------------------------------------------------------------*/

#ifndef directionalKSolidThermo_H
#define directionalKSolidThermo_H

#include "coordinateSystem.H"
#include "interpolatedSolidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class directionalKSolidThermo Declaration
\*---------------------------------------------------------------------------*/

class directionalKSolidThermo
:
    public interpolatedSolidThermo
{
    // Private data

        //- Thermal conductivity [W/(m.K)]
        volSymmTensorField directionalK_;

        //- Thermal conductivity vector
        Field<vector> KValues_;

        //- Coordinate system used for the directional properties
        coordinateSystem coordSys_;

        //- Transformation for cell centres
        volTensorField ccTransforms_;


    // Private Member Functions

        //- Transform principal values of symmTensor
        symmTensor transformPrincipal
        (
            const tensor& tt,
            const vector& st
         ) const;

        //- Transform principal values of symmTensor
        void transformField
        (
            symmTensorField& fld,
            const tensorField& tt,
            const vectorField& st
        ) const;

        //- Calculate properties
        void calculate();

        //- Initialize thermo
        void init();


public:

    //- Runtime type information
    TypeName("directionalKSolidThermo");


    // Constructors

        //- Construct from mesh
        directionalKSolidThermo(const fvMesh& mesh);

         //- Construct from mesh and dictionary
        directionalKSolidThermo(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~directionalKSolidThermo();


    // Member Functions

        //- Update properties
        virtual void correct();

        //- Access functions

            //- Thermal conductivity [W/m/K]
            virtual tmp<volScalarField> K() const
            {
                notImplemented("directionalKSolidThermo::K() const");
                return tmp<volScalarField>(NULL);
            }

            //- Thermal conductivity [W/m/K]
            virtual tmp<volSymmTensorField> directionalK() const;

         // Per patch calculation

            //- Thermal conductivity [W//m/K]
            virtual tmp<scalarField> K(const label patchI) const
            {
                notImplemented("directionalKSolidThermo::K(const label) const");
                return tmp<scalarField>(NULL);
            }

            //- Thermal conductivity [W//m/K]
            virtual tmp<symmTensorField> directionalK(const label) const;


        // I-O

            //- Write the directionalKSolidThermo properties
            virtual bool writeData(Ostream& os) const;

            //- Read the directionalKSolidThermo properties
            virtual bool read();

            //- Read the directionalKSolidThermo properties
            bool read(const dictionary& dict);

            //- Ostream Operator
            friend Ostream& operator<<
            (
                Ostream& os,
                const directionalKSolidThermo& s
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
