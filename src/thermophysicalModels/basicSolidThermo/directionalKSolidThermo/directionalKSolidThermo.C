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

#include "directionalKSolidThermo.H"
#include "addToRunTimeSelectionTable.H"
#include "transform.H"
#include "transformField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(directionalKSolidThermo, 0);
    addToRunTimeSelectionTable
    (
        basicSolidThermo,
        directionalKSolidThermo,
        mesh
    );

    addToRunTimeSelectionTable
    (
        basicSolidThermo,
        directionalKSolidThermo,
        dictionary
    );

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::directionalKSolidThermo::directionalKSolidThermo
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    interpolatedSolidThermo(mesh, typeName + "Coeffs", dict),
    directionalK_
    (
        IOobject
        (
            "K",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimEnergy/dimTime/(dimLength*dimTemperature)
    ),
    ccTransforms_
    (
        IOobject
        (
            "ccTransforms",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimLength
    )
{
    init();
}


Foam::directionalKSolidThermo::directionalKSolidThermo(const fvMesh& mesh)
:
    interpolatedSolidThermo(mesh, typeName + "Coeffs"),
    directionalK_
    (
        IOobject
        (
            "K",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimEnergy/dimTime/(dimLength*dimTemperature)
    ),
    ccTransforms_
    (
        IOobject
        (
            "ccTransforms",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimLength
    )
{
    init();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::directionalKSolidThermo::~directionalKSolidThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::directionalKSolidThermo::init()
{
    KValues_ = Field<vector>(subDict(typeName + "Coeffs").lookup("KValues"));

    // Determine transforms for cell centres
    forAll(mesh_.C(), cellI)
    {
        vector dir = mesh_.C()[cellI] - coordSys_.origin();
        dir /= mag(dir);

        // Define local coordinate system with
        // - e1 : axis from cc to centre
        // - e3 : rotation axis
        coordinateSystem cs
        (
            "cc",
            coordSys_.origin(),
            coordSys_.e3(),     //z',e3
            dir                 //x',e1
        );

        ccTransforms_[cellI] = cs.R();
    }

    forAll(mesh_.C().boundaryField(), patchI)
    {
        const fvPatchVectorField& patchC = mesh_.C().boundaryField()[patchI];
        fvPatchTensorField& patchT = ccTransforms_.boundaryField()[patchI];

        tensorField tc(patchT.size());
        forAll(tc, i)
        {
            vector dir = patchC[i] - coordSys_.origin();
            dir /= mag(dir);

            coordinateSystem cs
            (
                "cc",
                coordSys_.origin(),
                coordSys_.e3(),     //z',e3
                dir                 //x',e1
            );

            tc[i] = cs.R();
        }
        patchT = tc;
    }

    if (debug)
    {
        Info<< "directionalKSolidThermo : dumping converted Kxx, Kyy, Kzz"
            << endl;
        {
            volVectorField Kxx
            (
                IOobject
                (
                    "Kxx",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE,
                    false
                ),
                mesh_,
                dimless
            );
            Kxx.internalField() = transform
            (
                ccTransforms_.internalField(),
                vectorField
                (
                    ccTransforms_.internalField().size(),
                    point(1, 0, 0)
                )
            );
            forAll(Kxx.boundaryField(), patchI)
            {
                Kxx.boundaryField()[patchI] = transform
                (
                    ccTransforms_.boundaryField()[patchI],
                    vectorField
                    (
                        ccTransforms_.boundaryField()[patchI].size(),
                        point(1, 0, 0)
                    )
                );
            }
            Kxx.write();
        }
        {
            volVectorField Kyy
            (
                IOobject
                (
                    "Kyy",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE,
                    false
                ),
                mesh_,
                dimless
            );
            Kyy.internalField() = transform
            (
                ccTransforms_.internalField(),
                vectorField
                (
                    ccTransforms_.internalField().size(),
                    point(0, 1, 0)
                )
            );
            forAll(Kyy.boundaryField(), patchI)
            {
                Kyy.boundaryField()[patchI] = transform
                (
                    ccTransforms_.boundaryField()[patchI],
                    vectorField
                    (
                        ccTransforms_.boundaryField()[patchI].size(),
                        point(0, 1, 0)
                    )
                );
            }
            Kyy.write();
        }
        {
            volVectorField Kzz
            (
                IOobject
                (
                    "Kzz",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE,
                    false
                ),
                mesh_,
                dimless
            );
            Kzz.internalField() = transform
            (
                ccTransforms_.internalField(),
                vectorField
                (
                    ccTransforms_.internalField().size(),
                    point(0, 0, 1)
                )
            );
            forAll(Kzz.boundaryField(), patchI)
            {
                Kzz.boundaryField()[patchI] = transform
                (
                    ccTransforms_.boundaryField()[patchI],
                    vectorField
                    (
                        ccTransforms_.boundaryField()[patchI].size(),
                        point(0, 0, 1)
                    )
                );
            }
            Kzz.write();
        }
    }

    correct();
}


Foam::symmTensor Foam::directionalKSolidThermo::transformPrincipal
(
    const tensor& tt,
    const vector& st
) const
{
    return symmTensor
    (
        tt.xx()*st.x()*tt.xx()
      + tt.xy()*st.y()*tt.xy()
      + tt.xz()*st.z()*tt.xz(),

        tt.xx()*st.x()*tt.yx()
      + tt.xy()*st.y()*tt.yy()
      + tt.xz()*st.z()*tt.yz(),

        tt.xx()*st.x()*tt.zx()
      + tt.xy()*st.y()*tt.zy()
      + tt.xz()*st.z()*tt.zz(),

        tt.yx()*st.x()*tt.yx()
      + tt.yy()*st.y()*tt.yy()
      + tt.yz()*st.z()*tt.yz(),

        tt.yx()*st.x()*tt.zx()
      + tt.yy()*st.y()*tt.zy()
      + tt.yz()*st.z()*tt.zz(),

        tt.zx()*st.x()*tt.zx()
      + tt.zy()*st.y()*tt.zy()
      + tt.zz()*st.z()*tt.zz()
    );
}


void Foam::directionalKSolidThermo::transformField
(
    symmTensorField& fld,
    const tensorField& tt,
    const vectorField& st
) const
{
    fld.setSize(tt.size());
    forAll(fld, i)
    {
        fld[i] = transformPrincipal(tt[i], st[i]);
    }
}


void Foam::directionalKSolidThermo::correct()
{
    calculate();
    interpolatedSolidThermo::calculate();
}


Foam::tmp<Foam::volSymmTensorField>
Foam::directionalKSolidThermo::directionalK() const
{
    return directionalK_;
}


void Foam::directionalKSolidThermo::calculate()
{
    // Correct directionalK
    Field<vector> localK
    (
        interpolateXY
        (
            T_.internalField(),
            TValues_,
            KValues_
        )
    );

    // Transform into global coordinate system
    transformField
    (
        directionalK_.internalField(),
        ccTransforms_.internalField(),
        localK
    );

    forAll(directionalK_.boundaryField(), patchI)
    {
        directionalK_.boundaryField()[patchI] == this->directionalK(patchI)();
    }
}


Foam::tmp<Foam::symmTensorField> Foam::directionalKSolidThermo::directionalK
(
    const label patchI
) const
{
    const fvPatchScalarField& patchT = T_.boundaryField()[patchI];

    Field<vector> localK(interpolateXY(patchT, TValues_, KValues_));

    tmp<symmTensorField> tglobalK(new symmTensorField(localK.size()));
    transformField(tglobalK(), ccTransforms_.boundaryField()[patchI], localK);

    return tglobalK;
}


bool Foam::directionalKSolidThermo::read()
{
    return read(subDict(typeName + "Coeffs"));
}


bool Foam::directionalKSolidThermo::read(const dictionary& dict)
{
    coordSys_ = coordinateSystem(dict, mesh_);
    KValues_  = Field<vector>(subDict(typeName + "Coeffs").lookup("KValues"));
    return true;
}


bool Foam::directionalKSolidThermo::writeData(Ostream& os) const
{
    bool ok = interpolatedSolidThermo::writeData(os);
    os.writeKeyword("KValues") << KValues_ << token::END_STATEMENT << nl;
    return ok && os.good();
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const directionalKSolidThermo& s)
{
    s.writeData(os);
    return os;
}


// ************************************************************************* //
