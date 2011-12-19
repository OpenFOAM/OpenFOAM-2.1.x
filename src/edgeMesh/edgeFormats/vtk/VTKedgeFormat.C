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

#include "VTKedgeFormat.H"
#include "OFstream.H"
#include "clock.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fileFormats::VTKedgeFormat::writeHeader
(
    Ostream& os,
    const pointField& pointLst
)
{
    // Write header
    os  << "# vtk DataFile Version 2.0" << nl
        << "featureEdgeMesh written " << clock::dateTime().c_str() << nl
        << "ASCII" << nl
        << nl
        << "DATASET POLYDATA" << nl;

    // Write vertex coords
    os  << "POINTS " << pointLst.size() << " float" << nl;
    forAll(pointLst, ptI)
    {
        const point& pt = pointLst[ptI];

        os  << pt.x() << ' ' << pt.y() << ' ' << pt.z() << nl;
    }
}


void Foam::fileFormats::VTKedgeFormat::writeEdges
(
    Ostream& os,
    const UList<edge>& edgeLst
)
{
    os  << "LINES " << edgeLst.size() << ' ' << 3*edgeLst.size() << nl;

    forAll(edgeLst, edgeI)
    {
        const edge& e = edgeLst[edgeI];

        os  << "2 " << e[0] << ' ' << e[1] << nl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileFormats::VTKedgeFormat::VTKedgeFormat()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fileFormats::VTKedgeFormat::write
(
    const fileName& filename,
    const edgeMesh& eMesh
)
{
    OFstream os(filename);
    if (!os.good())
    {
        FatalErrorIn
        (
            "fileFormats::VTKedgeFormat::write"
            "(const fileName&, const edgeMesh&)"
        )
            << "Cannot open file for writing " << filename
            << exit(FatalError);
    }

    writeHeader(os, eMesh.points());
    writeEdges(os, eMesh.edges());
}


// ************************************************************************* //
