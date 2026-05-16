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
    ),
        dict_
    (
        mesh_.thisDb().lookupObject<IOdictionary>("fvSchemes")
    ),
        order_
    (
        // dict_.lookupOrAddDefault<word>("reconstructionOrder", "first")
        // dict_.lookup("reconstructionOrder")
        dict_.lookupOrDefault<word>("reconstructionOrder", "first")
    ),
    limiter_
    (
        // dict_.lookupOrAddDefault<word>("limiter", "noLimiter")
        // dict_.lookup("limiter")
        dict_.lookupOrDefault<word>("limiter", "no")
    )

{
    gradientSchemes::distanceMatrix(Ainv_);
    gradientSchemes::distanceMatrixLocal(AinvLocal_);

     if
    (
        order_ != "first" && order_ != "second"
    )
    {
        FatalErrorIn("gradientSchemes.H")
            << "Valid type entries are 'first' or 'second' "
            << "for reconstructionOrder"
            << abort(FatalError);
    }
     if
    (
        limiter_ != "no" && limiter_ != "yes"
    )
    {
        FatalErrorIn("gradientSchemes.H")
            << "Valid type entries are 'no' or 'yes' "
            << "for limiter"
            << abort(FatalError);
    }
    Info << " reconstruction scheme order: " << order_ <<endl;
    Info << " Limiter: " << limiter_ <<endl;
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
        forAll(U.boundaryField(), patchID)
        {

            // const fvPatch& curPatch = U.boundaryField()[patchi].patch();

            if (U.boundaryField()[patchID].coupled())
            {
                const fvPatch& p = mesh_.boundary()[patchID];
                // Better version of d-vectors: Zeljko Tukovic, 25/Apr/2010
                const vectorField pd = p.delta();

                // const vectorField X_nei
                // (
                //   X_.boundaryField()[patchID].patchNeighbourField()
                // );


                forAll(mesh_.boundary()[patchID], facei)
                {
                    const label& bCellID =
                        mesh_.boundaryMesh()[patchID].faceCells()[facei];

                    // const vector& d = X_nei[facei] - X_[bCellID];
                    const vector& d = pd[facei];
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
    const objectRegistry& db = mesh_.thisDb();
    const pointVectorField& lmN_ = db.lookupObject<pointVectorField> ("lmN");

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

            //! works only with 3D geometry. Further invistigation required
            // if (lmN_.boundaryField().types()[patchID] == "fixedValue")
            // {
            //     const label& faceID =
            //         mesh_.boundary()[patchID].patch().start() + facei;

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
    // #pragma message("Compiling OPENFOAM_NOT_EXTEND branch")

    Ainv.primitiveFieldRef() = inv(dCd.internalField());
#else
    // #pragma message("Compiling EXTEND branch")

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
                const fvPatch& curPatch = mesh_.boundary()[patchID];
                // distance between two cell centers accross coupled pathes
                const vectorField pd = curPatch.delta();
   

                const scalarField U_nei
                (
                  U.boundaryField()[patchID].patchNeighbourField()
                );

                forAll(mesh_.boundary()[patchID], facei)
                {
                    const label& bCellID =
                        mesh_.boundaryMesh()[patchID].faceCells()[facei];

                    const vector& d = pd[facei];

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
            //! This is changed after checking with Parallel run and it works
            // if (lmN_.boundaryField().types()[patchID] == "fixedValue")
            // {
            //     const label& faceID =
            //         mesh_.boundary()[patchID].patch().start() + facei;

            //     forAll(mesh_.faces()[faceID], nodei)
            //     {
            //         const label& nodeID = mesh_.faces()[faceID][nodei];
            //         vector d = mesh_.points()[nodeID] - mesh_.C()[bCellID];

            //         UgradX[bCellID] +=
            //             AinvLocal_[bCellID]
            //           & ((lmN_[nodeID].x() - U[bCellID].x())*d);

            //         UgradY[bCellID] +=
            //             AinvLocal_[bCellID]
            //           & ((lmN_[nodeID].y() - U[bCellID].y())*d);

            //         UgradZ[bCellID] +=
            //             AinvLocal_[bCellID]
            //           & ((lmN_[nodeID].z() - U[bCellID].z())*d);

            //         for (int i=0; i<7; i++)
            //         {
            //             d =
            //                 ((((i+1)*mesh_.points()[nodeID])
            //               + ((7-i)*mesh_.Cf().boundaryField()[patchID][facei]))/8.0)
            //               - mesh_.C()[bCellID];

            //             UgradX[bCellID] +=
            //                 AinvLocal_[bCellID]
            //               & ((lmN_[nodeID].x() - U[bCellID].x())*d);

            //             UgradY[bCellID] +=
            //                 AinvLocal_[bCellID]
            //               & ((lmN_[nodeID].y() - U[bCellID].y())*d);

            //             UgradZ[bCellID] +=
            //                 AinvLocal_[bCellID]
            //               & ((lmN_[nodeID].z() - U[bCellID].z())*d);
            //         }
            //     }
            // }
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
        
        else if (mesh_.boundary()[patchID].coupled())
        {

            const fvPatch& curPatch = mesh_.boundary()[patchID];
            // distance between two cell centers accross coupled pathes
            const vectorField pd = curPatch.delta();
            // distance between the patch face center and its owner cell center 
            vectorField pDeltaRLeft = curPatch.fvPatch::delta();
            // distance between the patch face center and its neighbore cell center 
            vectorField pDdeltaRRight = pDeltaRLeft - pd; 

            // tmp<vectorField> tmp_X_nei = X_.boundaryField()[patchID].patchNeighbourField();
            // const vectorField& X_nei = tmp_X_nei();

            tmp<scalarField> tmp_U_nei = U.boundaryField()[patchID].patchNeighbourField();
            const scalarField& U_nei = tmp_U_nei();
           
            tmp<vectorField> tmp_Ugrad_nei = Ugrad.boundaryField()[patchID].patchNeighbourField();
            const vectorField& Ugrad_nei = tmp_Ugrad_nei();

            forAll(mesh_.boundary()[patchID], facei)
            {
                const label& bCellID =
                    mesh_.boundaryMesh()[patchID].faceCells()[facei];

#ifdef OPENFOAM_NOT_EXTEND
                Up.boundaryFieldRef()[patchID][facei] =
                    U_nei[facei] + ( Ugrad_nei[facei] & pDdeltaRRight[facei]);

                Um.boundaryFieldRef()[patchID][facei] =
                    U[bCellID] + ( Ugrad[bCellID]& pDeltaRLeft[facei]);
#else
                Up.boundaryField()[patchID][facei] =
                    U_nei[facei] + ( Ugrad_nei[facei]& pDdeltaRRight[facei]);

                Um.boundaryField()[patchID][facei] =
                    U[bCellID] + ( Ugrad[bCellID]& pDeltaRLeft[facei]);
#endif


            }
        }
        else
        {
            forAll(mesh_.boundaryMesh()[patchID],facei)
            {
                const label& bCellID =
                    mesh_.boundaryMesh()[patchID].faceCells()[facei];

               //!this is a buge from original code
                // U.boundaryFieldRef()[patchID][facei] =
                //     U[bCellID] + ( Ugrad[bCellID]
                // & (XF_.boundaryField()[patchID][facei] - X_[bCellID]));


#ifdef OPENFOAM_NOT_EXTEND
                Um.boundaryFieldRef()[patchID][facei] =
                    U[bCellID] + ( Ugrad[bCellID]
                & (mesh_.Cf().boundaryField()[patchID][facei] - mesh_.C()[bCellID]));
#else
                Um.boundaryField()[patchID][facei] =
                    U[bCellID] + ( Ugrad[bCellID]
                & (mesh_.Cf().boundaryField()[patchID][facei] - mesh_.C()[bCellID]));
#endif
            
            
            
            }
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
        
        else if (mesh_.boundary()[patchID].coupled())
        {
            const fvPatch& curPatch = mesh_.boundary()[patchID];
            // distance between two cell centers accross coupled pathes
            const vectorField pd = curPatch.delta();
            // distance between the patch face center and its owner cell center 
            vectorField pDeltaRLeft = curPatch.fvPatch::delta();
            // distance between the patch face center and its neighbore cell center 
            vectorField pDdeltaRRight = pDeltaRLeft - pd; 
               
            // tmp<vectorField> tmp_X_nei = X_.boundaryField()[patchID].patchNeighbourField();
            // const vectorField& X_nei = tmp_X_nei();

            tmp<vectorField> tmp_U_nei = U.boundaryField()[patchID].patchNeighbourField();
            const vectorField& U_nei = tmp_U_nei();
           
            tmp<tensorField> tmp_Ugrad_nei = Ugrad.boundaryField()[patchID].patchNeighbourField();
            const tensorField& Ugrad_nei = tmp_Ugrad_nei();

            forAll(mesh_.boundary()[patchID], facei)
            {
                const label& bCellID =
                    mesh_.boundaryMesh()[patchID].faceCells()[facei];


#ifdef OPENFOAM_NOT_EXTEND
                Up.boundaryFieldRef()[patchID][facei] =
                    U_nei[facei] + ( Ugrad_nei[facei]& pDdeltaRRight[facei]);

                Um.boundaryFieldRef()[patchID][facei] =
                    U[bCellID] + ( Ugrad[bCellID]& pDeltaRLeft[facei]);
#else
                Up.boundaryField()[patchID][facei] =
                    U_nei[facei] + ( Ugrad_nei[facei]& pDdeltaRRight[facei]);

                Um.boundaryField()[patchID][facei] =
                    U[bCellID] + ( Ugrad[bCellID]& pDeltaRLeft[facei]);
#endif            
            
            }
        }
        else
        {
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
                    U[bCellID] + (Ugrad[bCellID]
                & (mesh_.Cf().boundaryField()[patchID][facei] - mesh_.C()[bCellID]));
#endif                
            
            }
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
        
 
        else if (mesh_.boundary()[patchID].coupled())
        {

                const fvPatch& curPatch = mesh_.boundary()[patchID];
                // distance between two cell centers accross coupled pathes
                const vectorField pd = curPatch.delta();
                // distance between the patch face center and its owner cell center 
                vectorField pDeltaRLeft = curPatch.fvPatch::delta();
                // distance between the patch face center and its neighbore cell center 
                vectorField pDdeltaRRight = pDeltaRLeft - pd;
                               
            // tmp<vectorField> tmp_X_nei = X_.boundaryField()[patchID].patchNeighbourField();
            // const vectorField& X_nei = tmp_X_nei();

            tmp<vectorField> tmp_Ux_nei = Ux.boundaryField()[patchID].patchNeighbourField();
            const vectorField& Ux_nei = tmp_Ux_nei();

            tmp<vectorField> tmp_Uy_nei = Uy.boundaryField()[patchID].patchNeighbourField();
            const vectorField& Uy_nei = tmp_Uy_nei();
            
            tmp<vectorField> tmp_Uz_nei = Uz.boundaryField()[patchID].patchNeighbourField();
            const vectorField& Uz_nei = tmp_Uz_nei();
           
            tmp<tensorField> tmp_Ugradx_nei = UxGrad.boundaryField()[patchID].patchNeighbourField();
            const tensorField& UxGrad_nei = tmp_Ugradx_nei();
           
            tmp<tensorField> tmp_Ugrady_nei = UyGrad.boundaryField()[patchID].patchNeighbourField();
            const tensorField& UyGrad_nei = tmp_Ugrady_nei();
           
            tmp<tensorField> tmp_Ugradz_nei = UzGrad.boundaryField()[patchID].patchNeighbourField();
            const tensorField& UzGrad_nei = tmp_Ugradz_nei();

            forAll(mesh_.boundary()[patchID], facei)
            {
                const label& bCellID =
                    mesh_.boundaryMesh()[patchID].faceCells()[facei];


                const vector& reconsX_nei =
                    Ux_nei[facei] + (UxGrad_nei[facei] & pDdeltaRRight[facei]);

                const vector& reconsY_nei =
                    Uy_nei[facei] + (UyGrad_nei[facei] & pDdeltaRRight[facei]);

                const vector& reconsZ_nei =
                    Uz_nei[facei] + (UzGrad_nei[facei] & pDdeltaRRight[facei]);


#ifdef OPENFOAM_NOT_EXTEND

                Up.boundaryFieldRef()[patchID][facei] =
                    tensor(reconsX_nei, reconsY_nei, reconsZ_nei);
#else

                Up.boundaryField()[patchID][facei] =
                    tensor(reconsX_nei, reconsY_nei, reconsZ_nei);
#endif  
                //-----------------------------------------------------
                const vector& reconsX =
                    Ux[bCellID] + (UxGrad[bCellID] & pDeltaRLeft[facei]);

                const vector& reconsY =
                    Uy[bCellID] + (UyGrad[bCellID] & pDeltaRLeft[facei]);

                const vector& reconsZ =
                    Uz[bCellID] + (UzGrad[bCellID] & pDeltaRLeft[facei]);

                // U.boundaryFieldRef()[patchID][facei] =
                //     tensor(reconsX, reconsY, reconsZ);

#ifdef OPENFOAM_NOT_EXTEND
                Um.boundaryFieldRef()[patchID][facei] =
                    tensor(reconsX, reconsY, reconsZ);
#else
                Um.boundaryField()[patchID][facei] =
                    tensor(reconsX, reconsY, reconsZ);
#endif  
                
            }
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

            // U.boundaryFieldRef()[patchID][facei] =
            //     tensor(reconsX, reconsY, reconsZ);

#ifdef OPENFOAM_NOT_EXTEND
            Um.boundaryFieldRef()[patchID][facei] =
                tensor(reconsX, reconsY, reconsZ);
#else
            Um.boundaryField()[patchID][facei] =
                tensor(reconsX, reconsY, reconsZ);
#endif 
        }
    }

    tvf.clear();
    tsf.clear();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//!-----------------------------------------------------------------------
//!-----------------------------------------------------------------------
//!-----------------------------------------------------------------------

inline void gradientSchemes::limiter
(
    scalar& lim,
    const scalar& cellVolume,
    const scalar& deltaOneMax,
    const scalar& deltaOneMin,
    const scalar& extrapolate
)
// {
//     if(limiterType == "BarthJesperson")
//     {
//         BarthJespersonLimiter
//         (
//             cellVolume,
//             deltaOneMax,
//             deltaOneMin,
//             extrapolate
//         )
//     }
//     else if(limiterType == "VekrantKrishnan")
//     (
//         VekrantKrishnan
//         (
//             cellVolume,
//             deltaOneMax,
//             deltaOneMin,
//             extrapolate
//         )
        
//     ) 
//     else if(limiterType == "VekrantKrishnan2")
//     (
//         VekrantKrishnan2
//         (
//             cellVolume,
//             deltaOneMax,
//             deltaOneMin,
//             extrapolate
//         )
        
//     ) 
// }
// {
//     if (mag(extrapolate) < SMALL)
//     {
//         return;
//     }

//     if (extrapolate - deltaOneMax > SMALL)
//     {
//         lim = min(lim, deltaOneMax/extrapolate);
//     }
//     else if (extrapolate - deltaOneMin < -SMALL)
//     {
//         lim = min(lim, deltaOneMin/extrapolate);
//     }
// }

// {
//     scalar k = 7;
//     scalar epsilonSquare = pow3(k)*cellVolume;

//     if (mag(extrapolate) < SMALL)
//     {
//         // Limiter remains unchanged
//         return;
//     }
//     else if (extrapolate > VSMALL)
//     {
//         lim = max
//         (
//             0,
//             min
//             (
//                 (
//                     (sqr(deltaOneMax) + epsilonSquare)*extrapolate
//                     + 2*sqr(extrapolate)*deltaOneMax
//                 )/
//                 stabilise
//                 (
//                     extrapolate*
//                     (
//                         sqr(deltaOneMax)
//                         + 2*sqr(extrapolate)
//                         + deltaOneMax*extrapolate
//                         + epsilonSquare
//                     ),
//                     SMALL
//                 ),
//                 lim
//             )
//         );
//     }
//     else if (extrapolate < VSMALL)
//     {
//         lim = max
//         (
//             0,
//             min
//             (
//                 (
//                     (sqr(deltaOneMin) + epsilonSquare)*extrapolate
//                     + 2*sqr(extrapolate)*deltaOneMin
//                 )/
//                 stabilise
//                 (
//                     extrapolate*
//                     (
//                         sqr(deltaOneMin)
//                         + 2*sqr(extrapolate)
//                         + deltaOneMin*extrapolate
//                         + epsilonSquare
//                     ),
//                     SMALL
//                 ),
//                 lim
//             )
//         );
//     }
// }
// {
//     scalar r;

//     if (extrapolate > SMALL)
//     {
//         r = deltaOneMax/extrapolate;
//     }
//     else if (extrapolate < -SMALL)
//     {
//         r = deltaOneMin/extrapolate;
//     }
//     else
//     {
//         return;
//     }

//     // scalar phi = (sqr(r) + r) / (sqr(r) + 1);
//     scalar phi =  max((r + 1)/(r + 1/stabilise(r, SMALL)), 0);

//     // lim = min(lim, max(0, phi));
//     lim =  max(0, phi);
// }
{

    if (mag(extrapolate) < SMALL)
    {
        // Limiter remains unchanged
        return;
    }
    else if  (extrapolate - deltaOneMax > SMALL)
    {

       scalar y= deltaOneMax/extrapolate;
       scalar phiY = (sqr(y) + 2*y)/stabilise( sqr(y) + y + 2,SMALL);
        lim = max
        (
            0,
            min
            (
                phiY,
                lim
            )
        );
    }
    else if (extrapolate - deltaOneMin < -SMALL)
    {
       scalar y = deltaOneMin/extrapolate;

       scalar phiY = (sqr(y) + 2*y)/stabilise( sqr(y) + y + 2,SMALL);

        lim = max
        (
            0,
            min
            (
                phiY,
                lim
            )
        );
    }
}

// {
//      scalar k_ = 1.5;
//     if
//     (
//         deltaOneMax - deltaOneMin < SMALL
//         || mag(extrapolate) < SMALL
//     )
//     {
//         return;
//     }

//     scalar y;

//     if (extrapolate > 0)
//     {
//         y = deltaOneMax/extrapolate;
//     }
//     else
//     {
//         y = deltaOneMin/extrapolate;
//     }

//     if (y < k_)
//     {
//         const scalar C2 = (3 - 2*k_)/sqr(k_);
//         const scalar C3 = -1/(3*sqr(k_)) - 2/(3*k_)*C2;

//         lim = min(lim, y + C2*sqr(y) + C3*pow3(y));
//     }
// }

    // {
    //     scalar r = 1;

    //         if (mag(extrapolate) < SMALL)
    //         {
    //             return;
    //         }

    //         if (extrapolate > SMALL)
    //         {
    //             r =  deltaOneMax/extrapolate;
    //         }
    //         else if (extrapolate < VSMALL)
    //         {
    //             r =  deltaOneMin/extrapolate;
    //         }
 
    //     lim = (r + mag(r))/(1 + mag(r));

    //     // lim = min(lim, phi_r);
    // }

    // {
    //         scalar epsilonSquare =
    //             sqr(epsilonPrime_*(deltaOneMax - deltaOneMin));

    //         if (extrapolate > VSMALL)
    //         {
    //             lim = max
    //             (
    //                 0,
    //                 min
    //                 (
    //                     (
    //                         (sqr(deltaOneMax) + epsilonSquare)*extrapolate
    //                       + 2*sqr(extrapolate)*deltaOneMax
    //                     )/
    //                     stabilise
    //                     (
    //                         extrapolate*
    //                         (
    //                             sqr(deltaOneMax)
    //                           + 2*sqr(extrapolate)
    //                           + deltaOneMax*extrapolate
    //                           + epsilonSquare
    //                         ),
    //                         SMALL
    //                     ),
    //                     lim
    //                 )
    //             );
    //         }
    //         else if (extrapolate < VSMALL)
    //         {
    //             lim = max
    //             (
    //                 0,
    //                 min
    //                 (
    //                     (
    //                         (sqr(deltaOneMin) + epsilonSquare)*extrapolate
    //                       + 2*sqr(extrapolate)*deltaOneMin
    //                     )/
    //                     stabilise
    //                     (
    //                         extrapolate*
    //                         (
    //                             sqr(deltaOneMin)
    //                           + 2*sqr(extrapolate)
    //                           + deltaOneMin*extrapolate
    //                           + epsilonSquare
    //                         ),
    //                         SMALL
    //                     ),
    //                     lim
    //                 )
    //             );
    //         }
    //         else
    //         {
    //             // Limiter remains unchanged
    //         }
    //     }

void gradientSchemes::reconstruct
(
    GeometricField<scalar, fvPatchField, volMesh>& U,
    const GeometricField<vector, fvPatchField, volMesh>& Ugrad,
    GeometricField<scalar, fvsPatchField, surfaceMesh>& Um,
    GeometricField<scalar, fvsPatchField, surfaceMesh>& Up,
    GeometricField<scalar, fvPatchField, volMesh>& phi

)
{
    // ============================================================
    // First-order fallback
    // ============================================================
    if (order_ != "second")
    {
        forAll(mesh_.cells(), cellI)
        {
            phi[cellI] = 0.0;
        }
    }
    else if (limiter_ == "no")
    {
        forAll(mesh_.cells(), cellI)
        {
            phi[cellI] = 1.0;
        }
    }
    else if (limiter_ == "yes")
    {
        // ========================================================
        // Step 1: Compute neighbourhood extrema (Algorithm 3.1-1)
        // ========================================================
        // Collect min and max values from neighbourhood
        GeometricField<scalar, fvPatchField, volMesh> UMinValue
        (
            "UMinValue",
            U
        );

       GeometricField<scalar, fvPatchField, volMesh> UMaxValue
        (
            "UMaxValue",
            U
        );

        const labelUList& owner = mesh_.owner();
        const labelUList& neighbour = mesh_.neighbour();
        
        Field<scalar>& UMinIn = UMinValue.internalField();
        Field<scalar>& UMaxIn = UMaxValue.internalField();
        
        forAll(owner, faceI)
        {
            const label own = owner[faceI];
            const label nei = neighbour[faceI];

            UMinIn[own] = min(UMinIn[own], U[nei]);
            UMinIn[nei] = min(UMinIn[nei], U[own]);

            UMaxIn[own] = max(UMaxIn[own], U[nei]);
            UMaxIn[nei] = max(UMaxIn[nei], U[own]);
        }
                  // Coupled boundaries
        forAll (mesh_.boundary(), patchI)
        {
            if (mesh_.boundary()[patchI].coupled())
            {
                const scalarField UNei =
                    U.boundaryField()[patchI].patchNeighbourField();

                const labelList& fc =
                    U.boundaryField()[patchI].patch().faceCells();

                forAll (fc, faceI)
                {
                    const label& curFC = fc[faceI];

                    // min value
                    UMinIn[curFC] =
                        min(UMinIn[curFC], UNei[faceI]);

                    // max value
                    UMaxIn[curFC] =
                        max(UMaxIn[curFC], UNei[faceI]);
                }
            }
        }

        const DimensionedField<scalar, volMesh>& cellVolume = mesh_.V();
        const volVectorField& cellCentre = mesh_.C();
        const surfaceVectorField& faceCentre = mesh_.Cf();

        scalarField& phiLimiterIn = phi.internalField();
        const vectorField& gradPhiIn = Ugrad.internalField();

        // Compute limiter values, internal faces
        forAll (owner, faceI)
        {
            const label& own = owner[faceI];
            const label& nei = neighbour[faceI];

            vector deltaRLeft = faceCentre[faceI] - cellCentre[own];
            vector deltaRRight = faceCentre[faceI] - cellCentre[nei];

                // Find minimal limiter value in each cell

            // Owner side
            limiter
            (
                phiLimiterIn[own],
                cellVolume[own],
                UMaxIn[own] - phi[own],
                UMinIn[own] - phi[own],
                (deltaRLeft & gradPhiIn[own])
            );

            // Neighbour side
            limiter
            (
                phiLimiterIn[nei],
                cellVolume[nei],
                UMaxIn[nei] - phi[nei],
                UMinIn[nei] - phi[nei],
                (deltaRRight & gradPhiIn[nei])
            );
        }
        // Coupled boundaries
        forAll (U.boundaryField(), patchI)
        {
            if (U.boundaryField()[patchI].coupled())
            {
                // Get patch
                const fvPatch& p = U.boundaryField()[patchI].patch();

                const scalarField pNei =
                    U.boundaryField()[patchI].patchNeighbourField();

                const labelList& fc = p.faceCells();

                const vectorField deltaR = p.Cf() - p.Cn();

                // Get gradients
                const vectorField gradPhiLeft =
                    Ugrad.boundaryField()[patchI].patchInternalField();

                const vectorField gradPhiRight =
                    Ugrad.boundaryField()[patchI].patchNeighbourField();

                // Find minimal limiter value in each cell
                forAll (fc, faceI)
                {
                    const label& curFC = fc[faceI];

                    limiter
                    (
                        phiLimiterIn[curFC],
                        cellVolume[curFC],
                        UMaxIn[curFC] - phi[curFC],
                        UMinIn[curFC] - phi[curFC],
                        (deltaR[faceI] & gradPhiRight[faceI])
                    );
                }
            }
        }

        // Do parallel communication to correct limiter on
        // coupled boundaries
        phi.correctBoundaryConditions();

    }
    // Info << "max(phi)" <<max(phi.internalField())<<endl;
    // Info << "min(phi)" <<min(phi.internalField())<<endl;
    forAll(mesh_.owner(), faceID)
    {
        const label& ownID = mesh_.owner()[faceID];
        const label& neiID = mesh_.neighbour()[faceID];

        Um[faceID]  = U[ownID] + phi[ownID]* (Ugrad[ownID] & (mesh_.Cf()[faceID] - mesh_.C()[ownID]));
        Up[faceID]  = U[neiID] + phi[neiID]* (Ugrad[neiID] & (mesh_.Cf()[faceID] - mesh_.C()[neiID]));
    }

    forAll(mesh_.boundary(), patchID)
    {
        // Check if the boundary patch is of type "empty"
        if (mesh_.boundary()[patchID].type() == "empty")
        {
            continue;
        }
        
        else if (mesh_.boundary()[patchID].coupled())
        {

            const fvPatch& curPatch = mesh_.boundary()[patchID];
            // distance between two cell centers accross coupled pathes
            const vectorField pd = curPatch.delta();
            // distance between the patch face center and its owner cell center 
            vectorField pDeltaRLeft = curPatch.fvPatch::delta();
            // distance between the patch face center and its neighbore cell center 
            vectorField pDdeltaRRight = pDeltaRLeft - pd; 

            // tmp<vectorField> tmp_X_nei = X_.boundaryField()[patchID].patchNeighbourField();
            // const vectorField& X_nei = tmp_X_nei();

            tmp<scalarField> tmp_U_nei = U.boundaryField()[patchID].patchNeighbourField();
            const scalarField& U_nei = tmp_U_nei();

            tmp<scalarField> tmp_phi_nei = phi.boundaryField()[patchID].patchNeighbourField();
            const scalarField& phi_nei = tmp_phi_nei();
        
            tmp<vectorField> tmp_Ugrad_nei = Ugrad.boundaryField()[patchID].patchNeighbourField();
            const vectorField& Ugrad_nei = tmp_Ugrad_nei();

            forAll(mesh_.boundary()[patchID], facei)
            {
                const label& bCellID =
                    mesh_.boundaryMesh()[patchID].faceCells()[facei];

#ifdef OPENFOAM_NOT_EXTEND
                Up.boundaryFieldRef()[patchID][facei] =
                    U_nei[facei] + phi_nei[facei]*( Ugrad_nei[facei] & pDdeltaRRight[facei]);

                Um.boundaryFieldRef()[patchID][facei] =
                    U[bCellID] + phi[bCellID]*( Ugrad[bCellID]& pDeltaRLeft[facei]);
#else
                Up.boundaryField()[patchID][facei] =
                    U_nei[facei] + phi_nei[facei]*( Ugrad_nei[facei]& pDdeltaRRight[facei]);

                Um.boundaryField()[patchID][facei] =
                    U[bCellID] + phi[bCellID]*( Ugrad[bCellID]& pDeltaRLeft[facei]);
#endif


            }
        }
        else
        {
            forAll(mesh_.boundaryMesh()[patchID],facei)
            {
                const label& bCellID =
                    mesh_.boundaryMesh()[patchID].faceCells()[facei];

            //!this is a buge from original code
                // U.boundaryFieldRef()[patchID][facei] =
                //     U[bCellID] + ( Ugrad[bCellID]
                // & (XF_.boundaryField()[patchID][facei] - X_[bCellID]));


#ifdef OPENFOAM_NOT_EXTEND
                Um.boundaryFieldRef()[patchID][facei] =
                    U[bCellID] + phi[bCellID]*( Ugrad[bCellID]
                & (mesh_.Cf().boundaryField()[patchID][facei] - mesh_.C()[bCellID]));
                
                U.boundaryFieldRef()[patchI][facei]  =
                   U[bCellID] + phi[bCellID]*( Ugrad[bCellID]
                & (mesh_.Cf().boundaryField()[patchID][facei] - mesh_.C()[bCellID]));
                
#else
                Um.boundaryField()[patchID][facei] =
                    U[bCellID] + phi[bCellID]*( Ugrad[bCellID]
                & (mesh_.Cf().boundaryField()[patchID][facei] - mesh_.C()[bCellID]));
                
                U.boundaryField()[patchID][facei]  =
                   U[bCellID] + phi[bCellID]*( Ugrad[bCellID]
                & (mesh_.Cf().boundaryField()[patchID][facei] - mesh_.C()[bCellID]));
                
#endif
            
            
            
            }
        }
    }
    
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void gradientSchemes::reconstruct
(
    GeometricField<vector, fvPatchField, volMesh>& U,
    const GeometricField<tensor, fvPatchField, volMesh>& Ugrad,
    GeometricField<vector, fvsPatchField, surfaceMesh>& Um,
    GeometricField<vector, fvsPatchField, surfaceMesh>& Up,
    GeometricField<vector, fvPatchField, volMesh>& phi
)
{
  tmp<GeometricField<scalar, fvPatchField, volMesh> > tvf
    (
        new GeometricField<scalar, fvPatchField, volMesh>
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
    GeometricField<scalar, fvPatchField, volMesh> Ux = tvf();
    GeometricField<scalar, fvPatchField, volMesh> Uy = tvf();
    GeometricField<scalar, fvPatchField, volMesh> Uz = tvf();

    operations op(mesh_);
    Ux= U.component(0);
    Uy= U.component(1);
    Uz= U.component(2);

      tmp<GeometricField<scalar, fvPatchField, volMesh> > tvfPhi
    (
        new GeometricField<scalar, fvPatchField, volMesh>
        (
            IOobject
            (
                "reconstructPhi("+U.name()+')',
                U.instance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            phi.dimensions()
        )
    );

    GeometricField<scalar, fvPatchField, volMesh> phiX = tvfPhi();
    GeometricField<scalar, fvPatchField, volMesh> phiY = tvfPhi();
    GeometricField<scalar, fvPatchField, volMesh> phiZ = tvfPhi();

    phiX= phi.component(0);
    phiY= phi.component(1);
    phiZ= phi.component(2);




  tmp<GeometricField<vector, fvPatchField, volMesh> > tvf1
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
            Ugrad.dimensions()
        )
    );
    GeometricField<vector, fvPatchField, volMesh> UxGrad = tvf1();
    GeometricField<vector, fvPatchField, volMesh> UyGrad = tvf1();
    GeometricField<vector, fvPatchField, volMesh> UzGrad = tvf1();

    op.decomposeTensor(Ugrad, UxGrad, UyGrad, UzGrad);



    tmp<GeometricField<scalar, fvsPatchField, surfaceMesh> > tsf
    (
        new GeometricField<scalar, fvsPatchField, surfaceMesh>
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
    GeometricField<scalar, fvsPatchField, surfaceMesh> UmX = tsf();
    GeometricField<scalar, fvsPatchField, surfaceMesh> UmY = tsf();
    GeometricField<scalar, fvsPatchField, surfaceMesh> UmZ = tsf();
    GeometricField<scalar, fvsPatchField, surfaceMesh> UpX = tsf();
    GeometricField<scalar, fvsPatchField, surfaceMesh> UpY = tsf();
    GeometricField<scalar, fvsPatchField, surfaceMesh> UpZ = tsf();



    // Info << "------vector Reconstruct---------"<<endl;

    gradientSchemes::reconstruct(Ux, UxGrad, UmX, UpX, phiX);
    gradientSchemes::reconstruct(Uy, UyGrad, UmY, UpY, phiY);
    gradientSchemes::reconstruct(Uz, UzGrad, UmZ, UpZ, phiZ);
    // gradientSchemes::reconstruct(Ux, UxGrad, UmX, UpX);
    // gradientSchemes::reconstruct(Uy, UyGrad, UmY, UpY);
    // gradientSchemes::reconstruct(Uz, UzGrad, UmZ, UpZ);

    forAll(own_, faceID)
    {
        Um[faceID] = vector(UmX[faceID], UmY[faceID], UmZ[faceID]);
        Up[faceID] = vector(UpX[faceID], UpY[faceID], UpZ[faceID]);
    }

    forAll(mesh_.boundary(), patchID)
    {
        // Check if the boundary patch is of type "empty"
        if (mesh_.boundary()[patchID].type() == "empty")
        {
            continue;
        }
        
        else if (mesh_.boundary()[patchID].coupled())
        {
            const fvPatch& curPatch = mesh_.boundary()[patchID];
            // distance between two cell centers accross coupled pathes
            const vectorField pd = curPatch.delta();
            // distance between the patch face center and its owner cell center 
            vectorField pDeltaRLeft = curPatch.fvPatch::delta();
            // distance between the patch face center and its neighbore cell center 
            vectorField pDdeltaRRight = pDeltaRLeft - pd; 


            tmp<vectorField> tmp_U_nei = U.boundaryField()[patchID].patchNeighbourField();
            const vectorField& U_nei = tmp_U_nei();

            tmp<vectorField> tmp_phi_nei = phi.boundaryField()[patchID].patchNeighbourField();
            const vectorField& phi_nei = tmp_phi_nei();
           
            tmp<tensorField> tmp_Ugrad_nei = Ugrad.boundaryField()[patchID].patchNeighbourField();
            const tensorField& Ugrad_nei = tmp_Ugrad_nei();

            forAll(mesh_.boundary()[patchID], facei)
            {
                const label& bCellID =
                    mesh_.boundaryMesh()[patchID].faceCells()[facei];


#ifdef OPENFOAM_NOT_EXTEND
                Up.boundaryFieldRef()[patchID][facei] =
                     U_nei[facei] + cmptMultiply(phi_nei[facei], (Ugrad_nei[facei]& pDdeltaRRight[facei]));
                    

                Um.boundaryFieldRef()[patchID][facei] =
                    U[bCellID] + phi[bCellID]*( Ugrad[bCellID]& pDeltaRLeft[facei]);
#else
                Up.boundaryField()[patchID][facei] =
                    U[bCellID] + cmptMultiply(phi[bCellID], (Ugrad[bCellID] & pDeltaRLeft[facei]));
                    

                Um.boundaryField()[patchID][facei] =
                    U[bCellID] + cmptMultiply(phi[bCellID], (Ugrad[bCellID] & pDeltaRLeft[facei]));
#endif            
            
            }
        }

        else
        {
            forAll(mesh_.boundaryMesh()[patchID], facei)
            {
                const label& bCellID =
                    mesh_.boundaryMesh()[patchID].faceCells()[facei];
                
                const vector delta =
                        mesh_.Cf().boundaryField()[patchID][facei]
                        - mesh_.C()[bCellID];
                
                const scalar& reconsX =
                    Ux[bCellID] +  phiX[bCellID]*(UxGrad[bCellID] & delta);

                const scalar& reconsY =
                    Uy[bCellID] + phiY[bCellID]*(UyGrad[bCellID]&delta);

                const scalar& reconsZ =
                    Uz[bCellID] + phiZ[bCellID]*(UzGrad[bCellID]&delta);

    #ifdef OPENFOAM_NOT_EXTEND
                U.boundaryFieldRef()[patchID][facei] =
                    vector(reconsX, reconsY, reconsZ);

                Um.boundaryFieldRef()[patchID][facei] =
                    vector(reconsX, reconsY, reconsZ);
    #else
                U.boundaryField()[patchID][facei] =
                    vector(reconsX, reconsY, reconsZ);

                Um.boundaryField()[patchID][facei] =
                    vector(reconsX, reconsY, reconsZ);
    #endif 
            }
        }
    }

    tvf.clear();
    tsf.clear();
    tvf1.clear();
    tvfPhi.clear();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void gradientSchemes::reconstruct
(
    GeometricField<tensor, fvPatchField, volMesh>& U,
    const GeometricField<tensor, fvPatchField, volMesh>& UxGrad,
    const GeometricField<tensor, fvPatchField, volMesh>& UyGrad,
    const GeometricField<tensor, fvPatchField, volMesh>& UzGrad,
    GeometricField<tensor, fvsPatchField, surfaceMesh>& Um,
    GeometricField<tensor, fvsPatchField, surfaceMesh>& Up,
    GeometricField<tensor, fvPatchField, volMesh>& phi

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

        tmp<GeometricField<vector, fvPatchField, volMesh> > tvfphi
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
            phi.dimensions()
        )
    );

    GeometricField<vector, fvPatchField, volMesh> phiX = tvfphi();
    GeometricField<vector, fvPatchField, volMesh> phiY = tvfphi();
    GeometricField<vector, fvPatchField, volMesh> phiZ = tvfphi();
    op.decomposeTensor(phi, phiX, phiY, phiZ);

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

    // Info << "------Tendor Reconstruct---------"<<endl;
    gradientSchemes::reconstruct(Ux, UxGrad, UmX, UpX, phiX);
    gradientSchemes::reconstruct(Uy, UyGrad, UmY, UpY, phiY);
    gradientSchemes::reconstruct(Uz, UzGrad, UmZ, UpZ, phiZ);


    forAll(mesh_.cells(), cellID)
    {
        phi[cellID] = tensor(phiX[cellID], phiY[cellID], phiZ[cellID]);
    }

    forAll(own_, faceID)
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

        else if (mesh_.boundary()[patchID].coupled())
        {

            const fvPatch& curPatch = mesh_.boundary()[patchID];
            // distance between two cell centers accross coupled pathes
            const vectorField pd = curPatch.delta();
            // distance between the patch face center and its owner cell center 
            vectorField pDeltaRLeft = curPatch.fvPatch::delta();
            // distance between the patch face center and its neighbore cell center 
            vectorField pDdeltaRRight = pDeltaRLeft - pd;
                               
            // tmp<vectorField> tmp_X_nei = X_.boundaryField()[patchID].patchNeighbourField();
            // const vectorField& X_nei = tmp_X_nei();

            tmp<vectorField> tmp_Ux_nei = Ux.boundaryField()[patchID].patchNeighbourField();
            const vectorField& Ux_nei = tmp_Ux_nei();

            tmp<vectorField> tmp_Uy_nei = Uy.boundaryField()[patchID].patchNeighbourField();
            const vectorField& Uy_nei = tmp_Uy_nei();
            
            tmp<vectorField> tmp_Uz_nei = Uz.boundaryField()[patchID].patchNeighbourField();
            const vectorField& Uz_nei = tmp_Uz_nei();  

            tmp<vectorField> tmp_phiX_nei = phiX.boundaryField()[patchID].patchNeighbourField();
            const vectorField& phiX_nei = tmp_phiX_nei();
            
            tmp<vectorField> tmp_phiY_nei = phiY.boundaryField()[patchID].patchNeighbourField();
            const vectorField& phiY_nei = tmp_phiY_nei();
            
            tmp<vectorField> tmp_phiZ_nei = phiZ.boundaryField()[patchID].patchNeighbourField();
            const vectorField& phiZ_nei = tmp_phiZ_nei();
           
            tmp<tensorField> tmp_Ugradx_nei = UxGrad.boundaryField()[patchID].patchNeighbourField();
            const tensorField& UxGrad_nei = tmp_Ugradx_nei();
           
            tmp<tensorField> tmp_Ugrady_nei = UyGrad.boundaryField()[patchID].patchNeighbourField();
            const tensorField& UyGrad_nei = tmp_Ugrady_nei();
           
            tmp<tensorField> tmp_Ugradz_nei = UzGrad.boundaryField()[patchID].patchNeighbourField();
            const tensorField& UzGrad_nei = tmp_Ugradz_nei();

            forAll(mesh_.boundary()[patchID], facei)
            {
                const label& bCellID =
                    mesh_.boundaryMesh()[patchID].faceCells()[facei];

                const vector& reconsX_nei =
                    Ux_nei[facei]
                + cmptMultiply(phiX_nei[facei],
                    (UxGrad_nei[facei] & pDdeltaRRight[facei]));

                const vector& reconsY_nei =
                    Uy_nei[facei]
                + cmptMultiply(phiY_nei[facei],
                    (UyGrad_nei[facei] & pDdeltaRRight[facei]));

                const vector& reconsZ_nei =
                    Uz_nei[facei]
                + cmptMultiply(phiZ_nei[facei],
                    (UzGrad_nei[facei] & pDdeltaRRight[facei]));

                #ifdef OPENFOAM_NOT_EXTEND
                Up.boundaryFieldRef()[patchID][facei] =
                    tensor(reconsX_nei, reconsY_nei, reconsZ_nei);
                #else
                Up.boundaryField()[patchID][facei] =
                    tensor(reconsX_nei, reconsY_nei, reconsZ_nei);
                #endif

                //-----------------------------------------------------

                const vector& reconsX =
                    Ux[bCellID]
                + cmptMultiply(phiX[bCellID],
                    (UxGrad[bCellID] & pDeltaRLeft[facei]));

                const vector& reconsY =
                    Uy[bCellID]
                + cmptMultiply(phiY[bCellID],
                    (UyGrad[bCellID] & pDeltaRLeft[facei]));

                const vector& reconsZ =
                    Uz[bCellID]
                + cmptMultiply(phiZ[bCellID],
                    (UzGrad[bCellID] & pDeltaRLeft[facei]));

                #ifdef OPENFOAM_NOT_EXTEND
                Um.boundaryFieldRef()[patchID][facei] =
                    tensor(reconsX, reconsY, reconsZ);
                #else
                Um.boundaryField()[patchID][facei] =
                    tensor(reconsX, reconsY, reconsZ);
                #endif
            }
        }        
        else
        {
            forAll(mesh_.boundaryMesh()[patchID], facei)
            {
                const label& bCellID =
                    mesh_.boundaryMesh()[patchID].faceCells()[facei];

                const vector delta =
                    mesh_.Cf().boundaryField()[patchID][facei]
                    - mesh_.C()[bCellID];

                const vector& reconsX =
                    Ux[bCellID]
                + cmptMultiply(phiX[bCellID],
                    (UxGrad[bCellID] & delta));

                const vector& reconsY =
                    Uy[bCellID]
                + cmptMultiply(phiY[bCellID],
                    (UyGrad[bCellID] & delta));

                const vector& reconsZ =
                    Uz[bCellID]
                + cmptMultiply(phiZ[bCellID],
                    (UzGrad[bCellID] & delta));

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
    }

    tvf.clear();
    tsf.clear();
    tvfphi.clear();
}




} // End namespace Foam

// ************************************************************************* //
