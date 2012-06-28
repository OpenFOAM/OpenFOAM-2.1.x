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

\*---------------------------------------------------------------------------*/

#include "nuSgsURoughWallFunctionFvPatchScalarField.H"
#include "LESModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

scalar nuSgsURoughWallFunctionFvPatchScalarField::calcYPlusLam
(
    const scalar kappa,
    const scalar E
) const
{
    scalar ypl = 11.0;

    for (int i=0; i<10; i++)
    {
        ypl = log(E*ypl)/kappa;
    }

    return ypl;
}


tmp<scalarField> nuSgsURoughWallFunctionFvPatchScalarField::calcYPlus
(
    const scalarField& magUp
) const
{
    const label patchi = patch().index();

    const LESModel& lesModel = db().lookupObject<LESModel>("LESProperties");
    const tmp<volScalarField> tnu = lesModel.nu();
    const volScalarField& nu = tnu();
    const scalarField& nuw = nu.boundaryField()[patchi];

    const scalarField& ry =
        lesModel.U().mesh().nonOrthDeltaCoeffs().boundaryField()[patchi];

    tmp<scalarField> tyPlus(new scalarField(patch().size(), 0.0));
    scalarField& yPlus = tyPlus();

    if (roughnessHeight_ > 0.0)
    {
        // Rough Walls
        const scalar c_1 = 1/(90 - 2.25) + roughnessConstant_;
        static const scalar c_2 = 2.25/(90 - 2.25);
        static const scalar c_3 = 2.0*atan(1.0)/log(90/2.25);
        static const scalar c_4 = c_3*log(2.25);

        //if (KsPlusBasedOnYPlus_)
        {
            // If KsPlus is based on YPlus the extra term added to the law
            // of the wall will depend on yPlus
            forAll(yPlus, facei)
            {
                const scalar magUpara = magUp[facei];
                const scalar Re = magUpara/(ry[facei]*nuw[facei]);
                const scalar kappaRe = kappa_*Re;

                scalar yp = yPlusLam_;
                const scalar ryPlusLam = 1.0/yp;

                int iter = 0;
                scalar yPlusLast = 0.0;
                scalar dKsPlusdYPlus = roughnessHeight_*ry[facei];

                // Additional tuning parameter - nominally = 1
                dKsPlusdYPlus *= roughnessFactor_;

                do
                {
                    yPlusLast = yp;

                    // The non-dimensional roughness height
                    scalar KsPlus = yp*dKsPlusdYPlus;

                    // The extra term in the law-of-the-wall
                    scalar G = 0.0;

                    scalar yPlusGPrime = 0.0;

                    if (KsPlus >= 90)
                    {
                        const scalar t_1 = 1 + roughnessConstant_*KsPlus;
                        G = log(t_1);
                        yPlusGPrime = roughnessConstant_*KsPlus/t_1;
                    }
                    else if (KsPlus > 2.25)
                    {
                        const scalar t_1 = c_1*KsPlus - c_2;
                        const scalar t_2 = c_3*log(KsPlus) - c_4;
                        const scalar sint_2 = sin(t_2);
                        const scalar logt_1 = log(t_1);
                        G = logt_1*sint_2;
                        yPlusGPrime =
                            (c_1*sint_2*KsPlus/t_1) + (c_3*logt_1*cos(t_2));
                    }

                    scalar denom = 1.0 + log(E_*yp) - G - yPlusGPrime;
                    if (mag(denom) > VSMALL)
                    {
                        yp = (kappaRe + yp*(1 - yPlusGPrime))/denom;
                    }
                } while
                (
                    mag(ryPlusLam*(yp - yPlusLast)) > 0.0001
                 && ++iter < 10
                 && yp > VSMALL
                );

                yPlus[facei] = max(0.0, yp);
            }
        }
    }
    else
    {
        // Smooth Walls
        forAll(yPlus, facei)
        {
            const scalar magUpara = magUp[facei];
            const scalar Re = magUpara/(ry[facei]*nuw[facei]);
            const scalar kappaRe = kappa_*Re;

            scalar yp = yPlusLam_;
            const scalar ryPlusLam = 1.0/yp;

            int iter = 0;
            scalar yPlusLast = 0.0;

            do
            {
                yPlusLast = yp;
                yp = (kappaRe + yp)/(1.0 + log(E_*yp));

            } while (mag(ryPlusLam*(yp - yPlusLast)) > 0.0001 && ++iter < 10);

            yPlus[facei] = max(0.0, yp);
        }
    }

    return tyPlus;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nuSgsURoughWallFunctionFvPatchScalarField::
nuSgsURoughWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    kappa_(0.41),
    E_(9.8),
    yPlusLam_(calcYPlusLam(kappa_, E_)),
    roughnessHeight_(pTraits<scalar>::zero),
    roughnessConstant_(pTraits<scalar>::zero),
    roughnessFactor_(pTraits<scalar>::zero)
{}


nuSgsURoughWallFunctionFvPatchScalarField::
nuSgsURoughWallFunctionFvPatchScalarField
(
    const nuSgsURoughWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    kappa_(ptf.kappa_),
    E_(ptf.E_),
    yPlusLam_(calcYPlusLam(kappa_, E_)),
    roughnessHeight_(ptf.roughnessHeight_),
    roughnessConstant_(ptf.roughnessConstant_),
    roughnessFactor_(ptf.roughnessFactor_)
{}


nuSgsURoughWallFunctionFvPatchScalarField::
nuSgsURoughWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    E_(dict.lookupOrDefault<scalar>("E", 9.8)),
    yPlusLam_(calcYPlusLam(kappa_, E_)),
    roughnessHeight_(readScalar(dict.lookup("roughnessHeight"))),
    roughnessConstant_(readScalar(dict.lookup("roughnessConstant"))),
    roughnessFactor_(readScalar(dict.lookup("roughnessFactor")))
{}


nuSgsURoughWallFunctionFvPatchScalarField::
nuSgsURoughWallFunctionFvPatchScalarField
(
    const nuSgsURoughWallFunctionFvPatchScalarField& rwfpsf
)
:
    fixedValueFvPatchScalarField(rwfpsf),
    kappa_(rwfpsf.kappa_),
    E_(rwfpsf.E_),
    yPlusLam_(rwfpsf.yPlusLam_),
    roughnessHeight_(rwfpsf.roughnessHeight_),
    roughnessConstant_(rwfpsf.roughnessConstant_),
    roughnessFactor_(rwfpsf.roughnessFactor_)
{}


nuSgsURoughWallFunctionFvPatchScalarField::
nuSgsURoughWallFunctionFvPatchScalarField
(
    const nuSgsURoughWallFunctionFvPatchScalarField& rwfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(rwfpsf, iF),
    kappa_(rwfpsf.kappa_),
    E_(rwfpsf.E_),
    yPlusLam_(rwfpsf.yPlusLam_),
    roughnessHeight_(rwfpsf.roughnessHeight_),
    roughnessConstant_(rwfpsf.roughnessConstant_),
    roughnessFactor_(rwfpsf.roughnessFactor_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void nuSgsURoughWallFunctionFvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{
    const label patchi = patch().index();

    const LESModel& lesModel = db().lookupObject<LESModel>("LESProperties");
    const fvPatchVectorField& Uw = lesModel.U().boundaryField()[patchi];
    const tmp<volScalarField> tnu = lesModel.nu();
    const volScalarField& nu = tnu();
    const scalarField& nuw = nu.boundaryField()[patchi];

    const scalarField& ry =
        lesModel.U().mesh().nonOrthDeltaCoeffs().boundaryField()[patchi];

    // The flow velocity at the adjacent cell centre
    const scalarField magUp(mag(Uw.patchInternalField() - Uw));

    tmp<scalarField> tyPlus = calcYPlus(magUp);
    scalarField& yPlus = tyPlus();

    scalarField& nuSgsw = *this;

    forAll(yPlus, facei)
    {
        if (yPlus[facei] > yPlusLam_)
        {
            const scalar Re = magUp[facei]/(ry[facei]*nuw[facei]) + ROOTVSMALL;
            nuSgsw[facei] = nuw[facei]*(sqr(yPlus[facei])/Re - 1);
        }
    }

    fixedValueFvPatchScalarField::evaluate();
}


void nuSgsURoughWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeKeyword("kappa") << kappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("E") << E_ << token::END_STATEMENT << nl;
    os.writeKeyword("roughnessHeight")
        << roughnessHeight_ << token::END_STATEMENT << nl;
    os.writeKeyword("roughnessConstant")
        << roughnessConstant_ << token::END_STATEMENT << nl;
    os.writeKeyword("roughnessFactor")
        << roughnessFactor_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    nuSgsURoughWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
