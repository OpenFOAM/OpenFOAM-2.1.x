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

#include "IDDESDelta.H"
#include "addToRunTimeSelectionTable.H"
#include "wallDistReflection.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(IDDESDelta, 0);
    addToRunTimeSelectionTable(LESdelta, IDDESDelta, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::IDDESDelta::calcDelta()
{
    label nD = mesh().nGeometricD();

    const volScalarField& hmax = hmax_();

    // initialise wallNorm
    wallDistReflection wallNorm(mesh());

    const volVectorField& n = wallNorm.n();

    tmp<volScalarField> faceToFacenMax
    (
        new volScalarField
        (
            IOobject
            (
                "faceToFaceMax",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("zrero", dimLength, 0.0)
        )
    );

    const cellList& cells = mesh().cells();

    forAll(cells,cellI)
    {
        scalar deltaMaxTmp = 0.0;
        const labelList& cFaces = mesh().cells()[cellI];
        const vector nCell = n[cellI];
        forAll(cFaces, cFaceI)
        {
            label faceI = cFaces[cFaceI];
            const point& faceCentreI = mesh().faceCentres()[faceI];
            forAll(cFaces, cFaceJ)
            {
                label faceJ = cFaces[cFaceJ];
                const point& faceCentreJ = mesh().faceCentres()[faceJ];
                scalar tmp = (faceCentreJ - faceCentreI) & nCell;
                if (tmp > deltaMaxTmp)
                {
                    deltaMaxTmp = tmp;
                }
            }
        }
        faceToFacenMax()[cellI] = deltaMaxTmp;
    }

    if (nD == 3)
    {
        delta_.internalField() =
            deltaCoeff_
           *min
            (
                max
                (
                    max(cw_*wallDist(mesh()).y(), cw_*hmax),
                    faceToFacenMax()
                ),
                hmax
            );
    }
    else if (nD == 2)
    {
        WarningIn("IDDESDelta::calcDelta()")
            << "Case is 2D, LES is not strictly applicable\n"
            << endl;

        delta_.internalField() =
            deltaCoeff_
           *min
            (
                max(max(cw_*wallDist(mesh()).y(), cw_*hmax), faceToFacenMax()),
                hmax
            );
    }
    else
    {
        FatalErrorIn("IDDESDelta::calcDelta()")
            << "Case is not 3D or 2D, LES is not strictly applicable"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IDDESDelta::IDDESDelta
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dd
)
:
    LESdelta(name, mesh),
    hmax_(LESdelta::New("hmax", mesh, dd.parent())),
    deltaCoeff_
    (
        readScalar(dd.subDict(type()+"Coeffs").lookup("deltaCoeff"))
    ),
    cw_(0.15)
{
    dd.subDict(type() + "Coeffs").readIfPresent("cw", cw_);
    calcDelta();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::IDDESDelta::read(const dictionary& dd)
{
    dd.subDict(type() + "Coeffs").lookup("deltaCoeff") >> deltaCoeff_;
    calcDelta();
}


void Foam::IDDESDelta::correct()
{
    if (mesh_.changing())
    {
        calcDelta();
        hmax_().correct();
    }
}


// ************************************************************************* //
