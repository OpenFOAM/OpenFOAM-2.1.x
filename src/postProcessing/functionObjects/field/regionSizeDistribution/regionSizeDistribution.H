/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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
    Foam::regionSizeDistribution

Description
    Droplet size distribution calculation.

    Looks up a void-fraction (alpha) field and splits the mesh into regions
    based on where the field is below the threshold value. These
    regions ("droplets") can now be analysed.

    Regions:
    - (debug) write regions as a volScalarField
    - (debug) print for all regions the sum of volume and alpha*volume
    - print the regions connected to a user-defined set of patches.
      (in spray calculation these form the liquid core)
    - print the regions with too large volume. These are the 'background'
      regions.

    Fields:
    - write volScalarField alpha_liquidCore : alpha with outside liquid core
                                              set to 0.
                           alpha_background : alpha with outside background
                                              set to 0.

    Histogram:
    - determine histogram of diameter (given minDiameter, maxDiameter, nBins)
    - write graph of number of droplets per bin
    - write graph of sum, average and deviation of droplet volume per bin
    - write graph of sum, average and deviation of user-defined fields. For
      volVectorFields these are those of the 3 components and the magnitude.

    Sample input:

    functions
    {
        regionSizeDistribution
        {
            type            regionSizeDistribution;

            outputControl   timeStep;
            outputInterval  1;

            // Field to determine regions from
            field           alpha;
            // Patches that provide the liquid core
            patches         (inlet);
            // Delimit alpha regions
            threshold       0.4;

            // Fields to sample (no need to include alpha)
            fields          (p U);

            // Number of bins for histogram
            nBins           100;
            // Max droplet diameter
            maxDiameter     0.5e-4;
            //// Min droplet diameter (default is 0)
            //minDiameter     0;

            // Writing format
            setFormat       gnuplot;
        }
    }



SourceFiles
    regionSizeDistribution.C

\*---------------------------------------------------------------------------*/

#ifndef regionSizeDistribution_H
#define regionSizeDistribution_H

#include "pointFieldFwd.H"
#include "writer.H"
#include "Map.H"
#include "volFieldsFwd.H"
#include "wordReList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class mapPolyMesh;
class regionSplit;
class polyMesh;

/*---------------------------------------------------------------------------*\
                           Class regionSizeDistribution Declaration
\*---------------------------------------------------------------------------*/

class regionSizeDistribution
{
    // Private data

        //- Name of this set of regionSizeDistribution objects
        word name_;

        const objectRegistry& obr_;

        //- on/off switch
        bool active_;

        //- Name of field
        word alphaName_;

        //- Patches to walk from
        wordReList patchNames_;

        //- Clip value
        scalar threshold_;

        //- Maximum droplet diameter
        scalar maxDiam_;

        //- Minimum droplet diameter
        scalar minDiam_;

        //- Mumber of bins
        label nBins_;

        //- Names of fields to sample on regions
        wordReList fields_;

        //- Output formatter to write
        autoPtr<writer<scalar> > formatterPtr_;


    // Private Member Functions

        template<class Type>
        Map<Type> regionSum(const regionSplit&, const Field<Type>&) const;

        //- Get data in order
        template<class Type>
        List<Type> extractData(const UList<label>& keys, const Map<Type>&)
        const;

        void writeGraph
        (
            const coordSet& coords,
            const word& valueName,
            const scalarField& values
        ) const;

        //- Write volfields with the parts of alpha which are not
        //  droplets (liquidCore, backGround)
        void writeAlphaFields
        (
            const regionSplit& regions,
            const Map<label>& keepRegions,
            const Map<scalar>& regionVolume,
            const volScalarField& alpha
        ) const;

        //- Mark all regions starting at patches
        Map<label> findPatchRegions(const polyMesh&, const regionSplit&) const;

        //- Helper: divide if denom != 0
        static tmp<scalarField> divide(const scalarField&, const scalarField&);

        //- Given per-region data calculate per-bin average/deviation and graph
        void writeGraphs
        (
            const word& fieldName,              // name of field
            const labelList& indices,           // index of bin for each region
            const scalarField& sortedField,     // per region field data
            const scalarField& binCount,        // per bin number of regions
            const coordSet& coords              // graph data for bins
        ) const;

        //- Given per-cell data calculate per-bin average/deviation and graph
        void writeGraphs
        (
            const word& fieldName,              // name of field
            const scalarField& cellField,       // per cell field data

            const regionSplit& regions,         // per cell the region(=droplet)
            const labelList& sortedRegions,     // valid regions in sorted order
            const scalarField& sortedNormalisation,

            const labelList& indices,           // index of bin for each region
            const scalarField& binCount,        // per bin number of regions
            const coordSet& coords              // graph data for bins
        ) const;

        //- Disallow default bitwise copy construct
        regionSizeDistribution(const regionSizeDistribution&);

        //- Disallow default bitwise assignment
        void operator=(const regionSizeDistribution&);


public:

    //- Runtime type information
    TypeName("regionSizeDistribution");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        regionSizeDistribution
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    // Destructor

        virtual ~regionSizeDistribution();


    // Member Functions

        //- Return name of the set of regionSizeDistribution
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the regionSizeDistribution data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Calculate the regionSizeDistribution and write
        virtual void write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const pointField&)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "regionSizeDistributionTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
