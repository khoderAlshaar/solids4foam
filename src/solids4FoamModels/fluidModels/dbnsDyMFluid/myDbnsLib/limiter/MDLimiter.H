/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    MDLimiter

Description
    Multidimensional limiter for density based solver

Author
    Aleksandar Jemcov
    Updated by Hrvoje Jasak

\*---------------------------------------------------------------------------*/

#ifndef MDLimiter_H
#define MDLimiter_H

// #include "volFieldsFwd.H"
// #include "surfaceFieldsFwd.H"
// #include "surfaceFieldsFwd.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "BarthJespersenLimiter.H"
#include "VenkatakrishnanLimiter.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class MDLimiter Declaration
\*---------------------------------------------------------------------------*/

template
<
    class Type,
    typename LimiterFunctionType
>
class MDLimiter
{
public:

    // Public typedefs

    typedef Field<Type> FieldType;
    typedef GeometricField<Type, fvPatchField, volMesh> GeoFieldType;

    typedef Field<typename outerProduct<vector, Type>::type> GradFieldType;
    typedef GeometricField
    <
        typename outerProduct<vector, Type>::type,
        fvPatchField,
        volMesh
    > GeoGradFieldType;


private:

    // Private data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Field
        const GeoFieldType& phi_;

        //- Grad field
        const GeoGradFieldType& gradPhi_;

        //- Limiter
        GeoFieldType phiLimiter_;


public:

    // Constructors

        //- Construct from field and gradient field
        MDLimiter
        (
            const GeoFieldType& phi,
            const GeoGradFieldType& gradPhi
        )
        :
            mesh_(phi.mesh()),
            phi_(phi),
            gradPhi_(gradPhi),
            phiLimiter_
            (
                IOobject
                (
                    "phiLimiter",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensioned<Type>("pLimiter", dimless, pTraits<Type>::one)
            )
        {
            // Collect min and max values from neighbourhood
            GeoFieldType phiMinValue
            (
                "phiMinValue",
                phi
            );

            GeoFieldType phiMaxValue
            (
                "phiMaxValue",
                phi
            );

            const auto& owner = mesh_.owner();
            const auto& neighbour = mesh_.neighbour();

            // Field<Type>& phiMinIn = phiMinValue.internalFieldRef();
            // Field<Type>& phiMaxIn = phiMaxValue.internalFieldRef();

            forAll (owner, faceI)
            {
                const label& own = owner[faceI];
                const label& nei = neighbour[faceI];

                // min values
                phiMinValue[own] = min(phiMinValue[own], phi_[nei]);
                phiMinValue[nei] = min(phiMinValue[nei], phi_[own]);

                // max values
                phiMaxValue[own] = max(phiMaxValue[own], phi_[nei]);
                phiMaxValue[nei] = max(phiMaxValue[nei], phi_[own]);
            }

            // Coupled boundaries
            forAll (phi.boundaryField(), patchI)
            {
                if (phi.boundaryField()[patchI].coupled())
                {
                    const Field<Type> pNei =
                        phi.boundaryField()[patchI].patchNeighbourField();

                    const labelList& fc =
                        phi.boundaryField()[patchI].patch().faceCells();

                    forAll (fc, faceI)
                    {
                        const label& curFC = fc[faceI];

                        // min value
                        phiMinValue[curFC] =
                            min(phiMinValue[curFC], pNei[faceI]);

                        // max value
                        phiMaxValue[curFC] =
                            max(phiMaxValue[curFC], pNei[faceI]);
                    }
                }
            }
            
            phiMinValue.correctBoundaryConditions();
            phiMaxValue.correctBoundaryConditions();
            // Calculate limiter

            // Get geometrical information

            const DimensionedField<scalar, volMesh>& cellVolume = mesh_.V();
            const volVectorField& cellCentre = mesh_.C();
            const surfaceVectorField& faceCentre = mesh_.Cf();

            // Create limiter function
            LimiterFunctionType limitFunction;

            // FieldType& phiLimiterIn = phiLimiter_.internalField();

            const GradFieldType& gradPhiIn = gradPhi.internalField();

            // Compute limiter values, internal faces
            forAll (owner, faceI)
            {
                const label& own = owner[faceI];
                const label& nei = neighbour[faceI];

                vector deltaRLeft = faceCentre[faceI] - cellCentre[own];
                vector deltaRRight = faceCentre[faceI] - cellCentre[nei];

                // Find minimal limiter value in each cell

                // Owner side
                limitFunction.limiter
                (
                    phiLimiter_[own],
                    cellVolume[own],
                    phiMaxValue[own] - phi_[own],
                    phiMinValue[own] - phi_[own],
                    (deltaRLeft & gradPhiIn[own])
                );

                // Neighbour side
                limitFunction.limiter
                (
                    phiLimiter_[nei],
                    cellVolume[nei],
                    phiMaxValue[nei] - phi_[nei],
                    phiMinValue[nei] - phi_[nei],
                    (deltaRRight & gradPhiIn[nei])
                );
            }

            // Coupled boundaries
            forAll (phi.boundaryField(), patchI)
            {
                if (phi.boundaryField()[patchI].coupled())
                {
                    // Get patch
                    const fvPatch& p = phi.boundaryField()[patchI].patch();

                    const FieldType pNei =
                        phi.boundaryField()[patchI].patchNeighbourField();

                    const labelList& fc = p.faceCells();

                    const vectorField deltaR = p.Cf() - p.Cn();

                    // Get gradients
                    const GradFieldType gradPhiLeft =
                        gradPhi.boundaryField()[patchI].patchInternalField();

                    const GradFieldType gradPhiRight =
                        gradPhi.boundaryField()[patchI].patchNeighbourField();

                    // Find minimal limiter value in each cell
                    forAll (fc, faceI)
                    {
                        const label& curFC = fc[faceI];

                        limitFunction.limiter
                        (
                            phiLimiter_[curFC],
                            cellVolume[curFC],
                            phiMaxValue[curFC] - phi_[curFC],
                            phiMinValue[curFC] - phi_[curFC],
                            (deltaR[faceI] & gradPhiLeft[faceI])
                        );
                    }
                }
            }

            // Do parallel communication to correct limiter on
            // coupled boundaries
            phiLimiter_.correctBoundaryConditions();
        }


    // Destructor - default


    // Member functions

        //- Return limiter
        const GeoFieldType& phiLimiter() const
        {
            return phiLimiter_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
