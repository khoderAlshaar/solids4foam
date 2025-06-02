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

#include "gradientSchemes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(gradientSchemes, 0);


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

gradientSchemes::gradientSchemes
(
    const fvMesh& vm
)
:
    mesh_(vm),
    own_(mesh_.owner()),
    nei_(mesh_.neighbour()),
    X_(mesh_.C()),
    XF_(mesh_.Cf()),
    XN_(mesh_.points()),

    Ainv_
    (
        IOobject("Ainv", mesh_),
        mesh_,
        dimensionedTensor("Ainv", dimensionSet(0,2,0,0,0,0,0), tensor::zero)
    ),

    AinvLocal_
    (
        IOobject("AinvLocal", mesh_),
        mesh_,
        dimensionedTensor
        (
            "AinvLocal",
            dimensionSet(0,2,0,0,0,0,0),
            tensor::zero
        )
    )
{
    gradientSchemes::distanceMatrix(Ainv_);
    gradientSchemes::distanceMatrixLocal(AinvLocal_);
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

gradientSchemes::~gradientSchemes()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void gradientSchemes::distanceMatrix
(
    GeometricField<tensor, fvPatchField, volMesh>& U
)
{
    forAll(mesh_.owner(), faceID)
    {
        const label& ownCellID = mesh_.owner()[faceID];
        const label& neiCellID = mesh_.neighbour()[faceID];
        const vector& dOwn = mesh_.C()[neiCellID] - mesh_.C()[ownCellID];
        const vector& dNei  = mesh_.C()[ownCellID] - mesh_.C()[neiCellID];

        U[ownCellID] += dOwn*dOwn;
        U[neiCellID] += dNei*dNei;
    }

    if (Pstream::parRun())
    {
        forAll(mesh_.boundary(), patchID)
        {
            if (mesh_.boundary()[patchID].coupled())
            {
                const vectorField X_nei
                (
                  mesh_.C().boundaryField()[patchID].patchNeighbourField()
                );

                forAll(mesh_.boundary()[patchID], facei)
                {
                    const label& bCellID =
                        mesh_.boundaryMesh()[patchID].faceCells()[facei];

                    const vector& d = X_nei[facei] - mesh_.C()[bCellID];
                    U[bCellID] += d*d;
                }
            }
        }

        U.correctBoundaryConditions();
    }

#ifdef OPENFOAM_NOT_EXTEND
    U.primitiveFieldRef() = inv(U.internalField());
#else
    U.internalField() = inv(U.internalField());
#endif

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void gradientSchemes::distanceMatrixLocal
(
    GeometricField<tensor, fvPatchField, volMesh>& Ainv
) const
{
    // const objectRegistry& db = mesh_.thisDb();
    // const pointVectorField& lmN_ = db.lookupObject<pointVectorField> ("lmN");

    tmp<GeometricField<tensor, fvPatchField, volMesh> > tvf
    (
        new GeometricField<tensor, fvPatchField, volMesh>
        (
            IOobject("distanceMatrixLocal", mesh_),
            mesh_,
            dimensioned<tensor>("0", Ainv.dimensions(), pTraits<tensor>::zero)
        )
    );
    GeometricField<tensor, fvPatchField, volMesh> dCd = tvf();

    forAll(mesh_.owner(), faceID)
    {
        const label& ownID = mesh_.owner()[faceID];
        const label& neiID = mesh_.neighbour()[faceID];
        const vector& dOwn = mesh_.Cf()[faceID] - mesh_.C()[ownID];
        const vector& dNei = mesh_.Cf()[faceID] - mesh_.C()[neiID];

        dCd[ownID] += dOwn*dOwn;
        dCd[neiID] += dNei*dNei;
    }

    forAll(mesh_.boundary(), patchID)
    {
        // Check if the boundary patch is of type "empty"
        if (mesh_.boundary()[patchID].type() == "empty")
        {
            continue;
        }

        forAll(mesh_.boundary()[patchID], facei)
        {
            const label& bCellID =
                mesh_.boundaryMesh()[patchID].faceCells()[facei];

            vector d = mesh_.Cf().boundaryField()[patchID][facei] - mesh_.C()[bCellID];
            dCd[bCellID] += d*d;

            // if (lmN_.boundaryField().types()[patchID] == "fixedValue")
            // {
            //     const label& faceID =
            //         mesh_.boundary()[patchID].start() + facei;

            //     forAll(mesh_.faces()[faceID], nodei)
            //     {
            //         const label& nodeID = mesh_.faces()[faceID][nodei];

            //         d = mesh_.points()[nodeID] - mesh_.C()[bCellID];
            //         dCd[bCellID] += d * d;

            //         for (int i=0; i<7; i++)
            //         {
            //             d =
            //                 ((((i+1)*mesh_.points()[nodeID])
            //               + ((7 - i)*mesh_.Cf().boundaryField()[patchID][facei]))/8.0)
            //               - mesh_.C()[bCellID];
            //             dCd[bCellID] += d * d;
            //         }
            //     }
            // }
        }
    }


#ifdef OPENFOAM_NOT_EXTEND
    Ainv.primitiveFieldRef() = inv(dCd.internalField());
#else
    Ainv.internalField() = inv(dCd.internalField());
#endif

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

volVectorField gradientSchemes::gradient
(
    const GeometricField<scalar, fvPatchField, volMesh>& U
)   const
{
    tmp<GeometricField<vector, fvPatchField, volMesh> > tvf
    (
        new GeometricField<vector, fvPatchField, volMesh>
        (
            IOobject
            (
                "gradient("+U.name()+')',
                mesh_
            ),
            mesh_,
            dimensioned<vector>
            (
                "0",
                U.dimensions()/dimLength,
                pTraits<vector>::zero
            )
        )
    );
    GeometricField<vector, fvPatchField, volMesh> Ugrad = tvf();

    forAll(mesh_.owner(), faceID)
    {
        const label& cellID = mesh_.owner()[faceID];
        const label& neiID = mesh_.neighbour()[faceID];
        const vector& dcell = mesh_.C()[neiID] - mesh_.C()[cellID];
        const vector& dnei = mesh_.C()[cellID] - mesh_.C()[neiID];

        Ugrad[cellID] += Ainv_[cellID] & (U[neiID] - U[cellID])*dcell;
        Ugrad[neiID] += Ainv_[neiID] & (U[cellID] - U[neiID])*dnei;
    }

    if (Pstream::parRun())
    {
        forAll(mesh_.boundary(), patchID)
        {
            if (mesh_.boundary()[patchID].coupled())
            {
                const vectorField X_nei
                (
                  mesh_.C().boundaryField()[patchID].patchNeighbourField()
                );

                const scalarField U_nei
                (
                  U.boundaryField()[patchID].patchNeighbourField()
                );

                forAll(mesh_.boundary()[patchID], facei)
                {
                    const label& bCellID =
                        mesh_.boundaryMesh()[patchID].faceCells()[facei];

                    const vector& d = X_nei[facei] - mesh_.C()[bCellID];

                    Ugrad[bCellID] +=
                        Ainv_[bCellID] & (U_nei[facei]-U[bCellID])*d;
                }
            }
        }

        Ugrad.correctBoundaryConditions();
    }

    tvf.clear();

    return Ugrad;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

volTensorField gradientSchemes::gradient
(
    const GeometricField<vector, fvPatchField, volMesh>& U
)   const
{
    tmp<GeometricField<vector, fvPatchField, volMesh> > tvf
    (
        new GeometricField<vector, fvPatchField, volMesh>
        (
            IOobject
            (
                "gradient("+U.name()+')',
                mesh_
            ),
            mesh_,
            dimensioned<vector>
            (
                "0",
                U.dimensions()/dimLength,
                pTraits<vector>::zero
            )
        )
    );
    GeometricField<vector, fvPatchField, volMesh> UgradX = tvf();
    GeometricField<vector, fvPatchField, volMesh> UgradY = tvf();
    GeometricField<vector, fvPatchField, volMesh> UgradZ = tvf();

    UgradX = gradientSchemes::gradient(U.component(0));
    UgradY = gradientSchemes::gradient(U.component(1));
    UgradZ = gradientSchemes::gradient(U.component(2));

    tmp<GeometricField<tensor, fvPatchField, volMesh> > ttf
    (
        new GeometricField<tensor, fvPatchField, volMesh>
        (
            IOobject
            (
                "gradient("+U.name()+')',
                mesh_
            ),
            mesh_,
            dimensioned<tensor>
            (
                "0",
                U.dimensions()/dimLength,
                pTraits<tensor>::zero
            )
        )
    );
    GeometricField<tensor, fvPatchField, volMesh> Ugrad = ttf();

    forAll(mesh_.cells(), cellID)
    {
        Ugrad[cellID] = tensor(UgradX[cellID], UgradY[cellID], UgradZ[cellID]);
    }

    if( Pstream::parRun() )
    {
        Ugrad.correctBoundaryConditions();
    }

    tvf.clear();
    ttf.clear();

    return Ugrad;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void gradientSchemes::gradient
(
    const GeometricField<tensor, fvPatchField, volMesh>& U,
    GeometricField<tensor, fvPatchField, volMesh>& UgradX,
    GeometricField<tensor, fvPatchField, volMesh>& UgradY,
    GeometricField<tensor, fvPatchField, volMesh>& UgradZ
)   const
{
    tmp<GeometricField<vector, fvPatchField, volMesh> > tvf
    (
        new GeometricField<vector, fvPatchField, volMesh>
        (
            IOobject
            (
                "gradient("+U.name()+')',
                U.instance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensioned<vector>("0", U.dimensions(), pTraits<vector>::zero)
        )
    );
    GeometricField<vector, fvPatchField, volMesh> Ux = tvf();
    GeometricField<vector, fvPatchField, volMesh> Uy = tvf();
    GeometricField<vector, fvPatchField, volMesh> Uz = tvf();

    operations op(mesh_);
    op.decomposeTensor(U, Ux, Uy, Uz);

    if (Pstream::parRun())
    {
        Ux.correctBoundaryConditions();
        Uy.correctBoundaryConditions();
        Uz.correctBoundaryConditions();
    }

    UgradX = gradientSchemes::gradient(Ux);
    UgradY = gradientSchemes::gradient(Uy);
    UgradZ = gradientSchemes::gradient(Uz);

    tvf.clear();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

volTensorField gradientSchemes::localGradient
(
    const GeometricField<vector, fvPatchField, volMesh>& U,
    const GeometricField<vector, fvsPatchField, surfaceMesh>& Unei
) const
{
    const objectRegistry& db = mesh_.thisDb();
    const pointVectorField& lmN_ = db.lookupObject<pointVectorField> ("lmN");

    tmp<GeometricField<vector, fvPatchField, volMesh> > tvf
    (
        new GeometricField<vector, fvPatchField, volMesh>
        (
            IOobject
            (
                "gradient("+U.name()+')',
                mesh_
            ),
            mesh_,
            dimensioned<vector>
            (
                "0",
                U.dimensions()/dimLength,
                pTraits<vector>::zero
            )
        )
    );
    GeometricField<vector, fvPatchField, volMesh> UgradX = tvf();
    GeometricField<vector, fvPatchField, volMesh> UgradY = tvf();
    GeometricField<vector, fvPatchField, volMesh> UgradZ = tvf();

    tmp<GeometricField<tensor, fvPatchField, volMesh> > tvft
    (
        new GeometricField<tensor, fvPatchField, volMesh>
        (
            IOobject
            (
                "gradient("+U.name()+')',
                mesh_
            ),
            mesh_,
            dimensioned<tensor>
            (
                "0",
                U.dimensions()/dimLength,
                pTraits<tensor>::zero
            )
        )
    );
    GeometricField<tensor, fvPatchField, volMesh> Ugrad = tvft();

    forAll(mesh_.owner(), faceID)
    {
        const label& ownID = mesh_.owner()[faceID];
        const label& neiID = mesh_.neighbour()[faceID];
        const vector& dOwn = mesh_.Cf()[faceID] - mesh_.C()[ownID];
        const vector& dNei = mesh_.Cf()[faceID] - mesh_.C()[neiID];

        UgradX[ownID] += AinvLocal_[ownID] & ((Unei[faceID].x()-U[ownID].x())*dOwn);
        UgradY[ownID] += AinvLocal_[ownID] & ((Unei[faceID].y()-U[ownID].y())*dOwn);
        UgradZ[ownID] += AinvLocal_[ownID] & ((Unei[faceID].z()-U[ownID].z())*dOwn);

        UgradX[neiID] += AinvLocal_[neiID] & ((Unei[faceID].x()-U[neiID].x())*dNei);
        UgradY[neiID] += AinvLocal_[neiID] & ((Unei[faceID].y()-U[neiID].y())*dNei);
        UgradZ[neiID] += AinvLocal_[neiID] & ((Unei[faceID].z()-U[neiID].z())*dNei);
    }

    forAll(mesh_.boundary(), patchID)
    {
        // Check if the boundary patch is of type "empty"
        if (mesh_.boundary()[patchID].type() == "empty")
        {
            continue;
        }

        forAll(mesh_.boundary()[patchID], facei)
        {
            const label& bCellID =
                mesh_.boundaryMesh()[patchID].faceCells()[facei];

            vector d = mesh_.Cf().boundaryField()[patchID][facei] - mesh_.C()[bCellID];

            UgradX[bCellID] +=
                AinvLocal_[bCellID]
              & ((Unei.boundaryField()[patchID][facei].x() - U[bCellID].x())*d);

            UgradY[bCellID] +=
                AinvLocal_[bCellID]
              & ((Unei.boundaryField()[patchID][facei].y() - U[bCellID].y())*d);

            UgradZ[bCellID] +=
                AinvLocal_[bCellID]
              & ((Unei.boundaryField()[patchID][facei].z() - U[bCellID].z())*d);

            if (lmN_.boundaryField().types()[patchID] == "fixedValue")
            {
                const label& faceID =
                    mesh_.boundary()[patchID].patch().start()  + facei;

                forAll(mesh_.faces()[faceID], nodei)
                {
                    const label& nodeID = mesh_.faces()[faceID][nodei];
                    vector d = mesh_.points()[nodeID] - mesh_.C()[bCellID];

                    UgradX[bCellID] +=
                        AinvLocal_[bCellID]
                      & ((lmN_[nodeID].x() - U[bCellID].x())*d);

                    UgradY[bCellID] +=
                        AinvLocal_[bCellID]
                      & ((lmN_[nodeID].y() - U[bCellID].y())*d);

                    UgradZ[bCellID] +=
                        AinvLocal_[bCellID]
                      & ((lmN_[nodeID].z() - U[bCellID].z())*d);

                    for (int i=0; i<7; i++)
                    {
                        d =
                            ((((i+1)*mesh_.points()[nodeID])
                          + ((7-i)*mesh_.Cf().boundaryField()[patchID][facei]))/8.0)
                          - mesh_.C()[bCellID];

                        UgradX[bCellID] +=
                            AinvLocal_[bCellID]
                          & ((lmN_[nodeID].x() - U[bCellID].x())*d);

                        UgradY[bCellID] +=
                            AinvLocal_[bCellID]
                          & ((lmN_[nodeID].y() - U[bCellID].y())*d);

                        UgradZ[bCellID] +=
                            AinvLocal_[bCellID]
                          & ((lmN_[nodeID].z() - U[bCellID].z())*d);
                    }
                }
            }
        }
    }

    forAll(mesh_.cells(), cellID)
    {
        Ugrad[cellID] = tensor(UgradX[cellID], UgradY[cellID], UgradZ[cellID]);
    }

    tvf.clear();
    tvft.clear();

    return Ugrad;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void gradientSchemes::reconstruct
(
    GeometricField<scalar, fvPatchField, volMesh>& U,
    const GeometricField<vector, fvPatchField, volMesh>& Ugrad,
    GeometricField<scalar, fvsPatchField, surfaceMesh>& Um,
    GeometricField<scalar, fvsPatchField, surfaceMesh>& Up
)
{
    forAll(mesh_.owner(), faceID)
    {
        const label& ownID = mesh_.owner()[faceID];
        const label& neiID = mesh_.neighbour()[faceID];

        Um[faceID] = U[ownID] + (Ugrad[ownID] & (mesh_.Cf()[faceID] - mesh_.C()[ownID]));
        Up[faceID]  = U[neiID] + (Ugrad[neiID] & (mesh_.Cf()[faceID] - mesh_.C()[neiID]));
    }

    forAll(mesh_.boundary(), patchID)
    {
        // Check if the boundary patch is of type "empty"
        if (mesh_.boundary()[patchID].type() == "empty")
        {
            continue;
        }

        forAll(mesh_.boundaryMesh()[patchID],facei)
        {
            const label& bCellID =
                mesh_.boundaryMesh()[patchID].faceCells()[facei];

#ifdef OPENFOAM_NOT_EXTEND
            U.boundaryFieldRef()[patchID][facei] =
                U[bCellID] + ( Ugrad[bCellID]
              & (mesh_.Cf().boundaryField()[patchID][facei] - mesh_.C()[bCellID]));

            Um.boundaryFieldRef()[patchID][facei] =
                U[bCellID] + ( Ugrad[bCellID]
              & (mesh_.Cf().boundaryField()[patchID][facei] - mesh_.C()[bCellID]));
#else
            U.boundaryField()[patchID][facei] =
                U[bCellID]+ ( Ugrad[bCellID]
             & ( mesh_.Cf().boundaryField()[patchID][facei] - mesh_.C()[bCellID]));

            Um.boundaryField()[patchID][facei] =
                U[bCellID]+ ( Ugrad[bCellID]
             & ( mesh_.Cf().boundaryField()[patchID][facei] - mesh_.C()[bCellID]));
#endif
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void gradientSchemes::reconstruct
(
    GeometricField<vector, fvPatchField, volMesh>& U,
    const GeometricField<tensor, fvPatchField, volMesh>& Ugrad,
    GeometricField<vector, fvsPatchField, surfaceMesh>& Um,
    GeometricField<vector, fvsPatchField, surfaceMesh>& Up
)
{
    forAll(mesh_.owner(), faceID)
    {
        const label& ownID = mesh_.owner()[faceID];
        const label& neiID = mesh_.neighbour()[faceID];

        Um[faceID] = U[ownID] + (Ugrad[ownID] & (mesh_.Cf()[faceID] - mesh_.C()[ownID]));
        Up[faceID] = U[neiID] + (Ugrad[neiID] & (mesh_.Cf()[faceID] - mesh_.C()[neiID]));
    }

    forAll(mesh_.boundary(), patchID)
    {
        // Check if the boundary patch is of type "empty"
        if (mesh_.boundary()[patchID].type() == "empty")
        {
            continue;
        }

        forAll(mesh_.boundaryMesh()[patchID], facei)
        {
            const label& bCellID =
                mesh_.boundaryMesh()[patchID].faceCells()[facei];

#ifdef OPENFOAM_NOT_EXTEND
            Um.boundaryFieldRef()[patchID][facei] =
                U[bCellID] + (Ugrad[bCellID]
              & (mesh_.Cf().boundaryField()[patchID][facei] - mesh_.C()[bCellID]));
#else
            Um.boundaryField()[patchID][facei] =
                U[bCellID]+ ( Ugrad[bCellID]
             & ( mesh_.Cf().boundaryField()[patchID][facei] - mesh_.C()[bCellID]));
#endif
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void gradientSchemes::reconstruct
(
    GeometricField<tensor, fvPatchField, volMesh>& U,
    const GeometricField<tensor, fvPatchField, volMesh>& UxGrad,
    const GeometricField<tensor, fvPatchField, volMesh>& UyGrad,
    const GeometricField<tensor, fvPatchField, volMesh>& UzGrad,
    GeometricField<tensor, fvsPatchField, surfaceMesh>& Um,
    GeometricField<tensor, fvsPatchField, surfaceMesh>& Up
)
{
    tmp<GeometricField<vector, fvPatchField, volMesh> > tvf
    (
        new GeometricField<vector, fvPatchField, volMesh>
        (
            IOobject
            (
                "reconstruct("+U.name()+')',
                U.instance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            U.dimensions()
        )
    );
    GeometricField<vector, fvPatchField, volMesh> Ux = tvf();
    GeometricField<vector, fvPatchField, volMesh> Uy = tvf();
    GeometricField<vector, fvPatchField, volMesh> Uz = tvf();

    operations op(mesh_);
    op.decomposeTensor(U, Ux, Uy, Uz);

    tmp<GeometricField<vector, fvsPatchField, surfaceMesh> > tsf
    (
        new GeometricField<vector, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "reconstruct("+Um.name()+')',
                Um.instance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            Um.dimensions()
        )
    );
    GeometricField<vector, fvsPatchField, surfaceMesh> UmX = tsf();
    GeometricField<vector, fvsPatchField, surfaceMesh> UmY = tsf();
    GeometricField<vector, fvsPatchField, surfaceMesh> UmZ = tsf();
    GeometricField<vector, fvsPatchField, surfaceMesh> UpX = tsf();
    GeometricField<vector, fvsPatchField, surfaceMesh> UpY = tsf();
    GeometricField<vector, fvsPatchField, surfaceMesh> UpZ = tsf();

    op.decomposeTensor(Um, UmX, UmY, UmZ);
    op.decomposeTensor(Up, UpX, UpY, UpZ);

    gradientSchemes::reconstruct(Ux, UxGrad, UmX, UpX);
    gradientSchemes::reconstruct(Uy, UyGrad, UmY, UpY);
    gradientSchemes::reconstruct(Uz, UzGrad, UmZ, UpZ);

    forAll(mesh_.owner(), faceID)
    {
        Um[faceID] = tensor(UmX[faceID], UmY[faceID], UmZ[faceID]);
        Up[faceID] = tensor(UpX[faceID], UpY[faceID], UpZ[faceID]);
    }

    forAll(mesh_.boundary(), patchID)
    {
        // Check if the boundary patch is of type "empty"
        if (mesh_.boundary()[patchID].type() == "empty")
        {
            continue;
        }

        forAll(mesh_.boundaryMesh()[patchID], facei)
        {
            const label& bCellID =
                mesh_.boundaryMesh()[patchID].faceCells()[facei];

            const vector& reconsX =
                Ux[bCellID] + (UxGrad[bCellID]
              & (mesh_.Cf().boundaryField()[patchID][facei] - mesh_.C()[bCellID]));

            const vector& reconsY =
                Uy[bCellID] + (UyGrad[bCellID]
              & (mesh_.Cf().boundaryField()[patchID][facei] - mesh_.C()[bCellID]));

            const vector& reconsZ =
                Uz[bCellID] + (UzGrad[bCellID]
              & (mesh_.Cf().boundaryField()[patchID][facei] - mesh_.C()[bCellID]));

#ifdef OPENFOAM_NOT_EXTEND
            U.boundaryFieldRef()[patchID][facei] =
                tensor(reconsX, reconsY, reconsZ);

            Um.boundaryFieldRef()[patchID][facei] =
                tensor(reconsX, reconsY, reconsZ);
#else
            U.boundaryField()[patchID][facei] =
                tensor(reconsX, reconsY, reconsZ);

            Um.boundaryField()[patchID][facei] =
                tensor(reconsX, reconsY, reconsZ);
#endif
        }
    }

    tvf.clear();
    tsf.clear();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
