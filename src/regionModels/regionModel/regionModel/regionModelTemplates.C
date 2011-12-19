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

template<class Type>
void Foam::regionModels::regionModel::toPrimary
(
    const label regionPatchI,
    List<Type>& regionField
) const
{
    forAll(intCoupledPatchIDs_, i)
    {
        if (intCoupledPatchIDs_[i] == regionPatchI)
        {
            const mappedPatchBase& mpb =
                refCast<const mappedPatchBase>
                (
                    regionMesh().boundaryMesh()[regionPatchI]
                );
            mpb.reverseDistribute(regionField);
            return;
        }
    }

    FatalErrorIn("const void toPrimary(const label, List<Type>&) const")
        << "Region patch ID " << regionPatchI << " not found in region mesh"
        << abort(FatalError);
}


template<class Type>
void Foam::regionModels::regionModel::toRegion
(
    const label regionPatchI,
    List<Type>& primaryField
) const
{
    forAll(intCoupledPatchIDs_, i)
    {
        if (intCoupledPatchIDs_[i] == regionPatchI)
        {
            const mappedPatchBase& mpb =
                refCast<const mappedPatchBase>
                (
                    regionMesh().boundaryMesh()[regionPatchI]
                );
            mpb.distribute(primaryField);
            return;
        }
    }

    FatalErrorIn("const void toRegion(const label, List<Type>&) const")
        << "Region patch ID " << regionPatchI << " not found in region mesh"
        << abort(FatalError);
}


template<class Type, class BinaryOp>
void Foam::regionModels::regionModel::toPrimary
(
    const label regionPatchI,
    List<Type>& regionField,
    const BinaryOp& bop
) const
{
    forAll(intCoupledPatchIDs_, i)
    {
        if (intCoupledPatchIDs_[i] == regionPatchI)
        {
            const mappedPatchBase& mpb =
                refCast<const mappedPatchBase>
                (
                    regionMesh().boundaryMesh()[regionPatchI]
                );
            mpb.reverseDistribute(regionField, bop);
            return;
        }
    }

    FatalErrorIn
    (
        "const void toPrimary"
        "("
            "const label, "
            "List<Type>&, "
            "const BinaryOp&"
        ") const"
    )   << "Region patch ID " << regionPatchI << " not found in region mesh"
        << abort(FatalError);
}


template<class Type, class BinaryOp>
void Foam::regionModels::regionModel::toRegion
(
    const label regionPatchI,
    List<Type>& primaryField,
    const BinaryOp& bop
) const
{
    forAll(intCoupledPatchIDs_, i)
    {
        if (intCoupledPatchIDs_[i] == regionPatchI)
        {
            const mappedPatchBase& mpb =
                refCast<const mappedPatchBase>
                (
                    regionMesh().boundaryMesh()[regionPatchI]
                );
            mpb.distribute(primaryField, bop);
            return;
        }
    }

    FatalErrorIn
    (
        "const void toRegion(const label, List<Type>&, const BinaryOp&) const"
    )   << "Region patch ID " << regionPatchI << " not found in region mesh"
        << abort(FatalError);
}


// ************************************************************************* //
