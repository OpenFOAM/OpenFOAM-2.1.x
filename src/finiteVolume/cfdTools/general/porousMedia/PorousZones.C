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

#include "PorousZones.H"
#include "Time.H"
#include "volFields.H"
#include "fvm.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ZoneType>
template<class Type>
void Foam::PorousZones<ZoneType>::modifyDdt(fvMatrix<Type>& m) const
{
    forAll(*this, i)
    {
        this->operator[](i).modifyDdt(m);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ZoneType>
Foam::PorousZones<ZoneType>::PorousZones
(
    const fvMesh& mesh
)
:
    IOPtrList<ZoneType>
    (
        IOobject
        (
            "porousZones",
            mesh.time().constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        typename ZoneType::iNew(mesh)
    ),
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ZoneType>
template<class Type>
Foam::tmp<Foam::fvMatrix<Type> >
Foam::PorousZones<ZoneType>::ddt
(
    GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type> > tres = fvm::ddt(vf);
    modifyDdt(tres());
    return tres;
}


template<class ZoneType>
template<class Type>
Foam::tmp<Foam::fvMatrix<Type> >
Foam::PorousZones<ZoneType>::ddt
(
    const geometricOneField&,
    GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type> > tres = fvm::ddt(vf);
    modifyDdt(tres());
    return tres;
}


template<class ZoneType>
template<class Type>
Foam::tmp<Foam::fvMatrix<Type> >
Foam::PorousZones<ZoneType>::ddt
(
    const dimensionedScalar& rho,
    GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type> > tres = fvm::ddt(rho,vf);
    modifyDdt(tres());
    return tres;
}


template<class ZoneType>
template<class Type>
Foam::tmp<Foam::fvMatrix<Type> >
Foam::PorousZones<ZoneType>::ddt
(
    const volScalarField& rho,
    GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type> > tres = fvm::ddt(rho,vf);
    modifyDdt(tres());
    return tres;
}

template<class ZoneType>
void Foam::PorousZones<ZoneType>::addResistance(fvVectorMatrix& UEqn) const
{
    forAll(*this, i)
    {
        this->operator[](i).addResistance(UEqn);
    }
}


template<class ZoneType>
void Foam::PorousZones<ZoneType>::addResistance
(
    fvVectorMatrix& UEqn,
    const volScalarField& rho,
    const volScalarField& mu
) const
{
    forAll(*this, i)
    {
        this->operator[](i).addResistance(UEqn, rho, mu);
    }
}


template<class ZoneType>
void Foam::PorousZones<ZoneType>::addResistance
(
    const fvVectorMatrix& UEqn,
    volTensorField& AU
) const
{
    // addResistance for each zone, delaying the correction of the
    // processor BCs of AU
    forAll(*this, i)
    {
        this->operator[](i).addResistance(UEqn, AU, false);
    }

    // Correct the boundary conditions of the tensorial diagonal to ensure
    // processor bounaries are correctly handled when AU^-1 is interpolated
    // for the pressure equation.
    AU.correctBoundaryConditions();
}


template<class ZoneType>
bool Foam::PorousZones<ZoneType>::readData(Istream& is)
{
    this->clear();

    IOPtrList<ZoneType> newLst
    (
        IOobject
        (
            "porousZones",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false     // Don't re-register new zones with objectRegistry
        ),
        typename ZoneType::iNew(mesh_)
    );

    this->transfer(newLst);

    return is.good();
}


template<class ZoneType>
bool Foam::PorousZones<ZoneType>::writeData(Ostream& os, bool subDict) const
{
    // Write size of list
    os << nl << this->size();

    // Write beginning of contents
    os << nl << token::BEGIN_LIST;

    // Write list contents
    forAll(*this, i)
    {
        os << nl;
        this->operator[](i).writeDict(os, subDict);
    }

    // Write end of contents
    os << token::END_LIST << nl;

    // Check state of IOstream
    return os.good();
}


// ************************************************************************* //
