/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "movingDisplacementNodalLinearMomentumPointPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "pointPatchFields.H"
#include "pointBoundaryMesh.H"
#include "pointMesh.H"
#ifdef OPENFOAM_NOT_EXTEND
    #include "Time.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

movingDisplacementNodalLinearMomentumPointPatchVectorField::
movingDisplacementNodalLinearMomentumPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchVectorField(p, iF),
    rho_(0.0),
    uMax_(vector::zero),
    tEnd_(0.0)
{}


movingDisplacementNodalLinearMomentumPointPatchVectorField::
movingDisplacementNodalLinearMomentumPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchVectorField(p, iF),
    rho_(dict.lookup("density")),
    uMax_(dict.lookup("displacement")),
    tEnd_(readScalar(dict.lookup("endTime")))
{
    fixedValuePointPatchVectorField::operator=(vectorField("value", dict, p.size()));
    updateCoeffs();
}


movingDisplacementNodalLinearMomentumPointPatchVectorField::
movingDisplacementNodalLinearMomentumPointPatchVectorField
(
    const movingDisplacementNodalLinearMomentumPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const PointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchVectorField(p, iF),
#ifdef OPENFOAM_NOT_EXTEND
    rho_(mapper(ptf.rho_)),
    uMax_(mapper(ptf.uMax_)),
    tEnd_(mapper(ptf.tEnd_))
#else
     rho_(ptf.rho_),
     uMax_(ptf.uMax_),
     tEnd_(ptf.tEnd_)
#endif
{}

#ifndef OPENFOAM_ORG
movingDisplacementNodalLinearMomentumPointPatchVectorField::
movingDisplacementNodalLinearMomentumPointPatchVectorField
(
    const movingDisplacementNodalLinearMomentumPointPatchVectorField& ptf
)
:
    fixedValuePointPatchVectorField(ptf),
    rho_(ptf.rho_),
    uMax_(ptf.uMax_),
    tEnd_(ptf.tEnd_)
{}
#endif

movingDisplacementNodalLinearMomentumPointPatchVectorField::
movingDisplacementNodalLinearMomentumPointPatchVectorField
(
    const movingDisplacementNodalLinearMomentumPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchVectorField(ptf, iF),
    rho_(ptf.rho_),
    uMax_(ptf.uMax_),
    tEnd_(ptf.tEnd_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void movingDisplacementNodalLinearMomentumPointPatchVectorField::rmap
(
    const pointPatchField<vector>& ptf,
    const labelList& addr
)
{
    fixedValuePointPatchVectorField::rmap(ptf, addr);
}


void movingDisplacementNodalLinearMomentumPointPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    vectorField lm_N(this->patchInternalField());

    const scalar& t = this->db().time().value();
    const scalar A = 3*10;
    const scalar B = 4*15/tEnd_;
    const scalar C = 5*6/pow(tEnd_,2);

    lm_N = rho_.value()*(uMax_/pow(tEnd_,3))*(A - B*t + C*t*t)*pow(t,2);

    this->operator==(lm_N);
    fixedValuePointPatchVectorField::updateCoeffs();
}


void movingDisplacementNodalLinearMomentumPointPatchVectorField::write
(
    Ostream& os
) const
{
    pointPatchVectorField::write(os);
    os.writeKeyword("density") << rho_ << token::END_STATEMENT << nl;
    os.writeKeyword("displacement") << uMax_ << token::END_STATEMENT << nl;
    os.writeKeyword("endTime") << tEnd_ << token::END_STATEMENT << nl;
    //writeEntry(os, "value", *this);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    movingDisplacementNodalLinearMomentumPointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
