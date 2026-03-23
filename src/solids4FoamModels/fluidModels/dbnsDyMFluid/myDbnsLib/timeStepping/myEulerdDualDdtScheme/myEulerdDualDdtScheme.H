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
    Foam::fv::myEulerdDualDdtScheme

Author
    Oliver Borm
    Aleksandar Jemcov

SourceFiles
    myEulerdDualDdtScheme.C

\*---------------------------------------------------------------------------*/

#ifndef myEulerdDualDdtScheme_H
#define myEulerdDualDdtScheme_H

#include "ddtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class myEulerdDualDdtScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class myEulerdDualDdtScheme
:
    public fv::ddtScheme<Type>
{
    // Private Member Functions

        //- Name of the physical time-step (field)
        word deltaTName_;

        //- Name of the pseudo time-step field
        word deltaTauName_;

        //- Name of old field
        word oldName_;
        //- current physical time-step
        scalar physDeltaT_;

        //- previous physical time-step
        scalar physDeltaT0_;

        //- Return the current time-step
        scalar deltaT_() const;

        //- Return the previous time-step
        scalar deltaT0_() const;

        //- Return the previous time-step or GREAT if the old timestep field
        //  wasn't available in which case Euler ddt is used
        template<class GeoField>
        scalar deltaT0_(const GeoField&) const;

        //- Disallow default bitwise copy construct
        myEulerdDualDdtScheme(const myEulerdDualDdtScheme&);

        //- Disallow default bitwise assignment
        void operator=(const myEulerdDualDdtScheme&);


public:

    //- Runtime type information
    TypeName("myEulerDual");


    // Constructors

        //- Construct from mesh
        myEulerdDualDdtScheme(const fvMesh& mesh)
        :
            ddtScheme<Type>(mesh),
            deltaTName_("physDeltaT"),
            deltaTauName_("CoDeltaT"),
            oldName_("_0"),
            physDeltaT_(0.1),
            physDeltaT0_(0.1)
        {}

        //- Construct from mesh and Istream
        myEulerdDualDdtScheme(const fvMesh& mesh, Istream& is)
        :
            ddtScheme<Type>(mesh, is),
            deltaTName_(is),
            deltaTauName_(is),
            oldName_("_0"),
            physDeltaT_(0.1),
            physDeltaT0_(0.1)
        {}


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return fv::ddtScheme<Type>::mesh();
        }

        tmp<GeometricField<Type, fvPatchField, volMesh> > fvcDdt
        (
            const dimensioned<Type>&
        );

        tmp<GeometricField<Type, fvPatchField, volMesh> > fvcDdt
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<GeometricField<Type, fvPatchField, volMesh> > fvcDdt
        (
            const dimensionedScalar&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<GeometricField<Type, fvPatchField, volMesh> > fvcDdt
        (
            const volScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<fvMatrix<Type> > fvmDdt
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<fvMatrix<Type> > fvmDdt
        (
            const dimensionedScalar&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<fvMatrix<Type> > fvmDdt
        (
            const volScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        typedef typename ddtScheme<Type>::fluxFieldType fluxFieldType;

        tmp<fluxFieldType> fvcDdtPhiCorr
        (
            const volScalarField& rA,
            const GeometricField<Type, fvPatchField, volMesh>& U,
            const fluxFieldType& phi
        );

        tmp<fluxFieldType> fvcDdtPhiCorr
        (
            const volScalarField& rA,
            const volScalarField& rho,
            const GeometricField<Type, fvPatchField, volMesh>& U,
            const fluxFieldType& phi
        );


        // Member functions for the new time consistent formulation
        tmp<fluxFieldType> fvcDdtConsistentPhiCorr
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>& faceU,
            const GeometricField<Type, fvPatchField, volMesh>& U,
            const surfaceScalarField& rAUf
        );


        tmp<surfaceScalarField> meshPhi
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );
};


template<>
tmp<surfaceScalarField> myEulerdDualDdtScheme<scalar>::fvcDdtPhiCorr
(
    const volScalarField& rA,
    const volScalarField& U,
    const surfaceScalarField& phi
);


template<>
tmp<surfaceScalarField> myEulerdDualDdtScheme<scalar>::fvcDdtPhiCorr
(
    const volScalarField& rA,
    const volScalarField& rho,
    const volScalarField& U,
    const surfaceScalarField& phi
);


template<>
tmp<surfaceScalarField> myEulerdDualDdtScheme<scalar>::fvcDdtConsistentPhiCorr
(
    const surfaceScalarField& faceU,
    const volScalarField& U,
    const surfaceScalarField& rAUf
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "myEulerdDualDdtScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
