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

#include "basicSource.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(basicSource, 0);
    defineRunTimeSelectionTable(basicSource, dictionary);

    template<> const char* NamedEnum
    <
        basicSource::selectionModeType,
        4
        >::names[] =
    {
        "points",
        "cellSet",
        "cellZone",
        "all"
    };

    const NamedEnum<basicSource::selectionModeType, 4>
        basicSource::selectionModeTypeNames_;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::basicSource::setSelection(const dictionary& dict)
{
    switch (selectionMode_)
    {
        case smPoints:
        {
            dict.lookup("points") >> points_;
            break;
        }
        case smCellSet:
        {
            dict.lookup("cellSet") >> cellSetName_;
            break;
        }
        case smCellZone:
        {
            dict.lookup("cellZone") >> cellSetName_;
            break;
        }
        case smAll:
        {
            break;
        }
        default:
        {
            FatalErrorIn
            (
                "basicSource::setSelection(const dictionary&)"
            )   << "Unknown selectionMode "
                << selectionModeTypeNames_[selectionMode_]
                << ". Valid selectionMode types are" << selectionModeTypeNames_
                << exit(FatalError);
        }
    }
}


void Foam::basicSource::setCellSet()
{
    Info<< incrIndent << indent << "Source: " << name_ << endl;
    switch (selectionMode_)
    {
        case smPoints:
        {
            Info<< indent << "- selecting cells using points" << endl;

            labelHashSet selectedCells;

            forAll(points_, i)
            {
                label cellI = mesh_.findCell(points_[i]);
                if (cellI >= 0)
                {
                    selectedCells.insert(cellI);
                }

                label globalCellI = returnReduce(cellI, maxOp<label>());
                if (globalCellI < 0)
                {
                    WarningIn("basicSource::setCellIds()")
                        << "Unable to find owner cell for point " << points_[i]
                        << endl;

                }

            }

            cells_ = selectedCells.toc();

            break;
        }
        case smCellSet:
        {
            Info<< indent << "- selecting cells using cellSet "
                << cellSetName_ << endl;

            cellSet selectedCells(mesh_, cellSetName_);
            cells_ = selectedCells.toc();

            break;
        }
        case smCellZone:
        {
            Info<< indent << "- selecting cells using cellZone "
                << cellSetName_ << endl;
            label zoneID = mesh_.cellZones().findZoneID(cellSetName_);
            if (zoneID == -1)
            {
                FatalErrorIn("basicSource::setCellIds()")
                    << "Cannot find cellZone " << cellSetName_ << endl
                    << "Valid cellZones are " << mesh_.cellZones().names()
                    << exit(FatalError);
            }
            cells_ = mesh_.cellZones()[zoneID];

            break;
        }
        case smAll:
        {
            Info<< indent << "- selecting all cells" << endl;
            cells_ = identity(mesh_.nCells());

            break;
        }
        default:
        {
            FatalErrorIn("basicSource::setCellIds()")
                << "Unknown selectionMode "
                << selectionModeTypeNames_[selectionMode_]
                << ". Valid selectionMode types are" << selectionModeTypeNames_
                << exit(FatalError);
        }
    }

    // Set volume information
    V_ = 0.0;
    forAll(cells_, i)
    {
        V_ += mesh_.V()[cells_[i]];
    }
    reduce(V_, sumOp<scalar>());

    Info<< indent << "- selected "
        << returnReduce(cells_.size(), sumOp<label>())
        << " cell(s) with volume " << V_ << nl << decrIndent << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicSource::basicSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    name_(name),
    mesh_(mesh),
    dict_(dict),
    coeffs_(dict.subDict(modelType + "Coeffs")),
    active_(readBool(dict_.lookup("active"))),
    timeStart_(readScalar(dict_.lookup("timeStart"))),
    duration_(readScalar(dict_.lookup("duration"))),
    selectionMode_
    (
        selectionModeTypeNames_.read(dict_.lookup("selectionMode"))
    ),
    cellSetName_("none"),
    V_(0.0),
    fieldNames_(),
    applied_()
{
    setSelection(dict_);

    setCellSet();
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::basicSource> Foam::basicSource::New
(
    const word& name,
    const dictionary& coeffs,
    const fvMesh& mesh
)
{
    word modelType(coeffs.lookup("type"));

    Info<< "Selecting source model type " << modelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "basicSource::New"
            "(const name&, const dictionary&, const fvMesh&)"
        )   << "Unknown Model type " << modelType
            << nl << nl
            << "Valid model types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<basicSource>(cstrIter()(name, modelType, coeffs, mesh));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::basicSource::isActive()
{
    if
    (
        active_
     && (mesh_.time().value() >= timeStart_)
     && (mesh_.time().value() <= timeEnd())
    )
    {
        // Update the cell set if the mesh is changing
        if (mesh_.changing())
        {
            setCellSet();
        }
        return true;
    }
    else
    {
        return false;
    }
}


Foam::label Foam::basicSource::applyToField(const word& fieldName) const
{
    forAll(fieldNames_, i)
    {
        if (fieldNames_[i] == fieldName)
        {
            return i;
        }
    }

    return -1;
}


void Foam::basicSource::checkApplied() const
{
    forAll(applied_, i)
    {
        if (!applied_[i])
        {
            WarningIn("void Foam::basicSource::checkApplied() const")
                << "Source " << name_ << " defined for field "
                << fieldNames_[i] << " but never used" << endl;
        }
    }
}


void Foam::basicSource::correct(volScalarField& fld)
{
    // do nothing
}


void Foam::basicSource::correct(volVectorField& fld)
{
    // do nothing
}


void Foam::basicSource::correct(volSphericalTensorField& fld)
{
    // do nothing
}


void Foam::basicSource::correct(volSymmTensorField& fld)
{
    // do nothing
}


void Foam::basicSource::correct(volTensorField& fld)
{
    // do nothing
}


void Foam::basicSource::addSup(fvMatrix<scalar>& eqn, const label fieldI)
{
    // do nothing
}


void Foam::basicSource::addSup(fvMatrix<vector>& eqn, const label fieldI)
{
    // do nothing
}


void Foam::basicSource::addSup
(
    fvMatrix<sphericalTensor>& eqn,
    const label fieldI
)
{
    // do nothing
}


void Foam::basicSource::addSup(fvMatrix<symmTensor>& eqn, const label fieldI)
{
    // do nothing
}


void Foam::basicSource::addSup(fvMatrix<tensor>& eqn, const label fieldI)
{
    // do nothing
}


void Foam::basicSource::setValue(fvMatrix<scalar>& eqn, const label fieldI)
{
    // do nothing
}


void Foam::basicSource::setValue(fvMatrix<vector>& eqn, const label fieldI)
{
    // do nothing
}


void Foam::basicSource::setValue
(
    fvMatrix<sphericalTensor>& eqn,
    const label fieldI
)
{
    // do nothing
}


void Foam::basicSource::setValue
(
    fvMatrix<symmTensor>& eqn,
    const label fieldI
)
{
    // do nothing
}


void Foam::basicSource::setValue(fvMatrix<tensor>& eqn, const label fieldI)
{
    // do nothing
}


// ************************************************************************* //
