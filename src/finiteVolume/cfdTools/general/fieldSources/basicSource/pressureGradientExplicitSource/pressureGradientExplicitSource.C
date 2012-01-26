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

#include "pressureGradientExplicitSource.H"
#include "fvMatrices.H"
#include "DimensionedField.H"
#include "IFstream.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pressureGradientExplicitSource, 0);

    addToRunTimeSelectionTable
    (
        basicSource,
        pressureGradientExplicitSource,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::pressureGradientExplicitSource::writeGradP() const
{
    // Only write on output time
    if (mesh_.time().outputTime())
    {
        IOdictionary propsDict
        (
            IOobject
            (
                name_ + "Properties",
                mesh_.time().timeName(),
                "uniform",
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            )
        );
        propsDict.add("gradient", gradP_);
        propsDict.regIOobject::write();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pressureGradientExplicitSource::pressureGradientExplicitSource
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    basicSource(sourceName, modelType, dict, mesh),
    Ubar_(coeffs_.lookup("Ubar")),
    gradPini_(coeffs_.lookup("gradPini")),
    gradP_(gradPini_),
    flowDir_(Ubar_/mag(Ubar_))
{
    coeffs_.lookup("fieldNames") >> fieldNames_;

    if (fieldNames_.size() != 1)
    {
        FatalErrorIn
        (
            "Foam::pressureGradientExplicitSource::"
            "pressureGradientExplicitSource"
            "("
                "onst word&, "
                "const word&, "
                "const dictionary&, "
                "const fvMesh&"
            ")"
        )   << "Source can only be applied to a single field.  Current "
            << "settings are:" << fieldNames_ << exit(FatalError);
    }

    applied_.setSize(fieldNames_.size(), false);

    // Read the initial pressure gradient from file if it exists
    IFstream propsFile
    (
        mesh_.time().timeName()/"uniform"/(name_ + "Properties")
    );

    if (propsFile.good())
    {
        Info<< "    Reading pressure gradient from file" << endl;
        dictionary propsDict(dictionary::null, propsFile);
        propsDict.lookup("gradient") >> gradP_;
    }

    Info<< "    Initial pressure gradient = " << gradP_ << nl << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pressureGradientExplicitSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    DimensionedField<vector, volMesh> Su
    (
        IOobject
        (
            name_ + fieldNames_[fieldI] + "Sup",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("zero", gradP_.dimensions(), vector::zero)
    );

    UIndirectList<vector>(Su, cells_) = flowDir_*gradP_.value();

    eqn += Su;
}


void Foam::pressureGradientExplicitSource::correct(volVectorField& U)
{
    const volScalarField& rAU =
        mesh_.lookupObject<volScalarField>("(1|A(" + U.name() + "))");

    // Integrate flow variables over cell set
    scalar magUbarAve = 0.0;
    scalar rAUave = 0.0;
    const scalarField& cv = mesh_.V();
    forAll(cells_, i)
    {
        label cellI = cells_[i];
        scalar volCell = cv[cellI];
        magUbarAve += (flowDir_ & U[cellI])*volCell;
        rAUave += rAU[cellI]*volCell;
    }

    // Collect across all processors
    reduce(magUbarAve, sumOp<scalar>());
    reduce(rAUave, sumOp<scalar>());

    // Volume averages
    magUbarAve /= V_;
    rAUave /= V_;

    // Calculate the pressure gradient increment needed to adjust the average
    // flow-rate to the desired value
    scalar gradPplus = (mag(Ubar_) - magUbarAve)/rAUave;

    // Apply correction to velocity field
    forAll(cells_, i)
    {
        label cellI = cells_[i];
        U[cellI] += flowDir_*rAU[cellI]*gradPplus;
    }

    // Update pressure gradient
    gradP_.value() += gradPplus;

    Info<< "Pressure gradient source: uncorrected Ubar = " << magUbarAve
        << ", pressure gradient = " << gradP_.value() << endl;

    writeGradP();
}


// ************************************************************************* //
