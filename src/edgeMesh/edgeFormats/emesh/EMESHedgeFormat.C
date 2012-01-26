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

#include "EMESHedgeFormat.H"
#include "IOobject.H"
#include "IFstream.H"
#include "clock.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileFormats::EMESHedgeFormat::EMESHedgeFormat
(
    const fileName& filename
)
{
    read(filename);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fileFormats::EMESHedgeFormat::read
(
    const fileName& filename
)
{
    clear();

    IFstream is(filename);
    if (!is.good())
    {
        FatalErrorIn
        (
            "fileFormats::EMESHedgeFormat::read(const fileName&)"
        )
            << "Cannot read file " << filename
            << exit(FatalError);
    }

    return read(is, this->storedPoints(), this->storedEdges());
}


bool Foam::fileFormats::EMESHedgeFormat::read
(
    Istream& is,
    pointField& pointLst,
    edgeList& edgeLst
)
{
    if (!is.good())
    {
        FatalErrorIn
        (
            "fileFormats::EMESHedgeFormat::read"
            "(Istream&, pointField&, edgeList&)"
        )
            << "read error "
            << exit(FatalError);
    }

    token firstToken(is);

    // swallow IOobject header
    if (!is.good())
    {
        FatalIOErrorIn
        (
            "fileFormats::EMESHedgeFormat::read"
            "(Istream&, pointField&, edgeList&)",
            is
        )
            << "First token could not be read" << nl
            << exit(FatalIOError);

        return false;
    }
    else if (firstToken.isWord() && firstToken.wordToken() == "FoamFile")
    {
        // read and discard
        dictionary headerDict(is);
    }
    else
    {
        is.putBack(firstToken);
    }

    // read points:
    is  >> pointLst;

    // read edges:
    is  >> edgeLst;

    return true;
}


Foam::Ostream& Foam::fileFormats::EMESHedgeFormat::write
(
    Ostream& os,
    const pointField& pointLst,
    const edgeList& edgeLst
)
{
    if (!os.good())
    {
        FatalErrorIn
        (
            "fileFormats::EMESHedgeFormat::write"
            "(Ostream&, const fileName&, const edgeMesh&)"
        )
            << "bad output stream " << os.name()
            << exit(FatalError);
    }

    os  << "\n// points:" << nl << pointLst << nl
        << "\n// edges:" << nl << edgeLst << nl;

    IOobject::writeDivider(os);

    // Check state of Ostream
    os.check
    (
        "EMESHedgeFormat::write"
        "(Ostream&, const pointField&, const edgeList&)"
    );

    return os;
}


void Foam::fileFormats::EMESHedgeFormat::write
(
    const fileName& filename,
    const edgeMesh& mesh
)
{
    OFstream os(filename);
    if (!os.good())
    {
        FatalErrorIn
        (
            "fileFormats::EMESHedgeFormat::write"
            "(const fileName&, const edgeMesh&)"
        )
            << "Cannot open file for writing " << filename
            << exit(FatalError);
    }


    // just emit some information until we get a nicer IOobject
    IOobject::writeBanner(os)
        << "FoamFile\n{\n"
        << "    version     " << os.version() << ";\n"
        << "    format      " << os.format() << ";\n"
        << "    class       " << "featureEdgeMesh" << ";\n"
        << "    note        " << "written " + clock::dateTime() << ";\n"
        << "    object      " << filename.name() << ";\n"
        << "}" << nl;

    IOobject::writeDivider(os);

    write(os, mesh.points(), mesh.edges());

    // Check state of Ostream
    os.check("EMESHedgeFormat::write(Ostream&)");
}


// ************************************************************************* //
