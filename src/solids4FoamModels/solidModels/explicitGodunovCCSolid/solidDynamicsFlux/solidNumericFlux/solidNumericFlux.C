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

\*---------------------------------------------------------------------------*/
  
#include "fvCFD.H"
#include "solidNumericFlux.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::solidNumericFlux::solidNumericFlux
(
    Time& runTime,
    const word& region,
    const dynamicFvMesh& mesh,
    volVectorField& lm,
    pointVectorField& lmN,
    volTensorField& F,
    volTensorField& P,
    solidMaterialModel& model,
    operations& op,
    mechanics& mech,
    gradientSchemes& grad
)
:

    mesh_(mesh),
    pFlux_( 
        Foam::solidFlux::New
        ( 
            mesh_, 
            mesh_.thisDb().lookupObject<IOdictionary>("fvSchemes")
        ) 
    ),
    lm_(lm),
    lmN_(lmN),
    F_(F),
    P_(P),

    model_(model),
    op_(op),
    mech_(mech),
    grad_(grad),
    // ----------------
    magSf_(mesh_.magSf()),
    Sf_(mesh_.Sf()),
    N_((Sf_ / mesh_.magSf())),
    n_(N_),

    rho_(model_.density()),

    Up_
    (
        IOobject("Up_", mesh_),
        mesh_,
        model.Up()
    ),

    Us_
    (
        IOobject("Us_", mesh_),
        mesh_,
        model.Us()
    ),
    lambda_
    (
        IOobject("lambda_", mesh_),
        mesh_,
        model.lambda()
    ),
    
    S_lm_(mech_.Smatrix_lm()),
    S_t_(mech_.Smatrix_t()),

    T1_
    (
            IOobject("T1_", mesh_),
            mesh_,
            dimensionedVector("T1_", dimensionSet(0,0,0,0,0,0,0), vector::zero)
    ),
    T2_
    (
        IOobject("T2_", mesh_),
        mesh_,
        dimensionedVector("T2_", dimensionSet(0,0,0,0,0,0,0), vector::zero)
    ),
    R_
    (
        IOobject("R", mesh_),
        mesh_,
        tensor::I
    ),
    RTranspos_
    (
        IOobject("RTranspos", mesh_),
        mesh_,
        tensor::I
    ),



// Reconstructed face fields (OWNED!)
    lm_P_
    (
        IOobject("lm_P", mesh_),
        mesh_,
        dimensionedVector("lm_P", lm_.dimensions(), vector::zero)

    ),
    lm_M_
    (
        IOobject("lm_M", mesh_),
        mesh_,
        dimensionedVector("lm_M", lm_.dimensions(), vector::zero)
    ),

    F_P_
    (
        IOobject("F_P", mesh_),
        mesh_,
        dimensionedTensor("F_P", F_.dimensions(), tensor::zero)
    ),
    F_M_
    (
        IOobject("F_M",  mesh_),
        mesh_,
        dimensionedTensor("F_M", F_.dimensions(), tensor::zero)

    ),

    P_P_
    (
        IOobject("P_P", mesh_),
        mesh_,
        dimensionedTensor("P_P", P_.dimensions(), tensor::zero)
    ),
    P_M_
    (
        IOobject("P_M", mesh_),
        mesh_,
        dimensionedTensor("P_M", P_.dimensions(), tensor::zero)

    ),

    // Tractions
    t_P_
    (
        IOobject("t_P", mesh_),
        mesh_,
        dimensionedVector("t_P", dimPressure, vector::zero)

    ),
    t_M_
    (
        IOobject("t_M",  mesh_),
        mesh_,
        dimensionedVector("t_M", dimPressure, vector::zero)
    ),
    lm_b_
    (
        IOobject
        (
            "lm_b",
             mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    t_b_
    (
        IOobject
        (
            "t_b",
             mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    lmFlux_
    (
        IOobject
        (
            "lmFlux",
             mesh_
        ),
        mesh_,
        dimensionedVector("lmFlux", dimensionSet(1,-1,-2,0,0,0,0), vector::zero)
    ),
    
    FFlux_
    (
        IOobject
        (
            "FFlux",
            mesh_
        ),
        mesh_,
        dimensionedTensor("FFlux", dimensionSet(0,1,-1,0,0,0,0), tensor::zero)
    )
{
    Info << "Hello from roeFlux constructor" << endl;
   setRotationalMatrix();


}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::solidNumericFlux::computeFlux()
{
    const auto& owner = mesh_.owner();
    // const auto& neighbour = mesh_.neighbour();

    // Get the face area vector
    // const surfaceVectorField& Sf = mesh_.Sf();
    // const surfaceScalarField& magSf = mesh_.magSf();

    //
    //!  Gradient evaluation
    // surfaceScalarField pos_(IOobject("pos", mesh_), mesh_, dimensionedScalar("one", dimless, 1.0));
    // surfaceScalarField neg_(IOobject("neg", mesh_), mesh_, dimensionedScalar("minusOne", dimless, -1.0));

    //  lm_P_ =  fvc::interpolate(lm_, pos_,"reconstruct(lm)") ;
    //  lm_M_ =  fvc::interpolate(lm_, neg_,"reconstruct(lm)") ;

    //  P_P_ =  fvc::interpolate(P_, pos_,"reconstruct(P)") ;
    //  P_M_ =  fvc::interpolate(P_, neg_,"reconstruct(P)") ;

    //  F_P_ =  fvc::interpolate(F_, pos_,"reconstruct(F)") ;
    //  F_M_ =  fvc::interpolate(F_, neg_,"reconstruct(F)") ;



    //decompoz of tensor F
    volVectorField Fx = op_.decomposeTensorX(F_);
    volVectorField Fy = op_.decomposeTensorY(F_);
    volVectorField Fz = op_.decomposeTensorZ(F_);

    // P = model.piola();
    volVectorField Px = op_.decomposeTensorX(P_);
    volVectorField Py = op_.decomposeTensorY(P_);
    volVectorField Pz = op_.decomposeTensorZ(P_);


    volTensorField lmGrad = grad_.gradient(lm_);

    volTensorField PxGrad = grad_.gradient(Px);
    volTensorField PyGrad = grad_.gradient(Py);
    volTensorField PzGrad = grad_.gradient(Pz);
    // compute F gradient
    volTensorField FxGrad = grad_.gradient(Fx);
    volTensorField FyGrad = grad_.gradient(Fy);
    volTensorField FzGrad = grad_.gradient(Fz);

    grad_.reconstruct(lm_, lmGrad, lm_M_, lm_P_);
    grad_.reconstruct(P_, PxGrad, PyGrad, PzGrad, P_M_, P_P_);
    grad_.reconstruct(F_, FxGrad, FyGrad, FzGrad, F_M_, F_P_);





//     // Riemann solver
    S_lm_ = mech_.Smatrix_lm();
    S_t_ = mech_.Smatrix_t();

   // Calculate fluxes at internal faces
    forAll (owner, faceI) 
    {
        // const label own = owner[faceI];
        // const label nei = neighbour[faceI];

        // t_M_[faceI] =   P_[own] & N_[faceI];

        // t_M_[faceI] =   P_M_[faceI] & N_[faceI];

        // calculate fluxes with reconstructed primitive variables at faces
        pFlux_ -> evaluateFlux
            (
                lmFlux_[faceI],
                FFlux_[faceI],
                lm_P_[faceI],
                lm_M_[faceI],
                F_P_[faceI],
                F_M_[faceI],
                P_P_[faceI],
                P_M_[faceI],
                R_[faceI],
                RTranspos_[faceI],
                rho_.value(),
                lambda_[faceI],
                Up_[faceI],
                Us_[faceI]
            );

            //     pFlux_ -> evaluateFlux
            // (
            //     lmFlux_[faceI],
            //     FFlux_[faceI],
            //     lm_[nei],
            //     lm_[own],
            //     F_[nei],
            //     F_[own],
            //     P_[nei],
            //     P_[own],
            //     R_[faceI],
            //     RTranspos_[faceI],
            //     rho_.value(),
            //     lambda_[faceI],
            //     Up_[faceI],
            //     Us_[faceI]
            // );
    }




    // Update boundary field and values
    forAll (lmFlux_.boundaryField(), patchi)
    {
        const fvPatch& curPatch = lm_.boundaryField()[patchi].patch();

        // Fluxes

#ifdef OPENFOAM_NOT_EXTEND
        fvsPatchVectorField& pLmFlux = lmFlux_.boundaryFieldRef()[patchi];
        fvsPatchTensorField& pFFlux = FFlux_.boundaryFieldRef()[patchi];
#else
        fvsPatchVectorField& pLmFlux = lmFlux_.boundaryField()[patchi];
        fvsPatchTensorField& pFFlux = FFlux_.boundaryField()[patchi];
#endif

        
    //    const labelList& fc =
    //                     lm_.boundaryField()[patchi].patch().faceCells();

        if (curPatch.coupled())
        {
            // Coupled patch
            // const vectorField plmLeft  =
            //     lm_.boundaryField()[patchi].patchInternalField();

            // const vectorField plmRight =
            //     lm_.boundaryField()[patchi].patchNeighbourField();

            // const tensorField pFLeft  =
            //     F_.boundaryField()[patchi].patchInternalField();

            // const tensorField pFRight =
            //     F_.boundaryField()[patchi].patchNeighbourField();

            // const tensorField pPLeft  =
            //     P_.boundaryField()[patchi].patchInternalField();

            // const tensorField pPRight =
            //     P_.boundaryField()[patchi].patchNeighbourField();



            // Patch fields
#ifdef OPENFOAM_NOT_EXTEND
            const fvsPatchVectorField& plm_P = lm_P_.boundaryFieldRef()[patchi];
            const fvsPatchTensorField& pP_P  = P_P_.boundaryFieldRef()[patchi];
            const fvsPatchTensorField& pF_P  = F_P_.boundaryFieldRef()[patchi];
            // const fvsPatchVectorField& pt_P  = t_P_.boundaryFieldRef()[patchi];
            
            const fvsPatchVectorField& plm_M = lm_M_.boundaryFieldRef()[patchi];
            const fvsPatchTensorField& pP_M = P_M_.boundaryFieldRef()[patchi];
            const fvsPatchTensorField& pF_M = F_M_.boundaryFieldRef()[patchi];
            // const fvsPatchVectorField& pt_M = t_M_.boundaryFieldRef()[patchi];

            const fvsPatchScalarField& pUp = Up_.boundaryFieldRef()[patchi];
            const fvsPatchScalarField& pUs = Us_.boundaryFieldRef()[patchi];
            const fvsPatchScalarField& pLmabda = lambda_.boundaryFieldRef()[patchi];
            
            // Face areas
            const fvsPatchTensorField& pR = R_.boundaryFieldRef()[patchi];
            const fvsPatchTensorField& pRT = RTranspos_.boundaryFieldRef()[patchi];

#else
            // Patch fields
            const fvsPatchVectorField& plm_P = lm_P_.boundaryField()[patchi];
            const fvsPatchTensorField& pP_P  = P_P_.boundaryField()[patchi];
            const fvsPatchTensorField& pF_P  = F_P_.boundaryField()[patchi];
            // const fvsPatchVectorField& pt_P  = t_P_.boundaryField()[patchi];
            
            const fvsPatchVectorField& plm_M = lm_M_.boundaryField()[patchi];
            const fvsPatchTensorField& pP_M = P_M_.boundaryField()[patchi];
            const fvsPatchTensorField& pF_M = F_M_.boundaryField()[patchi];
            // const fvsPatchVectorField& pt_M = t_M_.boundaryField()[patchi];

            const fvsPatchScalarField& pUp = Up_.boundaryField()[patchi];
            const fvsPatchScalarField& pUs = Us_.boundaryField()[patchi];
            const fvsPatchScalarField& pLmabda = lambda_.boundaryField()[patchi];
            
            // Face areas
            const fvsPatchTensorField& pR = R_.boundaryField()[patchi];
            const fvsPatchTensorField& pRT = RTranspos_.boundaryField()[patchi];
#endif
            forAll (curPatch, facei)
            {

                // const label& curFC = fc[facei];
                
            //    t_M_.boundaryField()[patchi][facei] =
            //                          P_[curFC] & N_.boundaryField()[patchi][facei];



// #ifdef OPENFOAM_NOT_EXTEND
//                t_M_.boundaryFieldRef()[patchi][facei] =
//                                       pPLeft[facei] & N_.boundaryField()[patchi][facei];
// #else
//                t_M_.boundaryField()[patchi][facei] =
//                                       pPLeft[facei] & N_.boundaryField()[patchi][facei];
// #endif

                    // const tensor& R = R_.boundaryField()[patchi][facei];
                    // const tensor& RT = RTranspos_.boundaryField()[patchi][facei];

                pFlux_ -> evaluateFlux
                (
                    pLmFlux[facei],
                    pFFlux[facei],
                    plm_P[facei],
                    plm_M[facei],
                    pF_P[facei],
                    pF_M[facei],
                    pP_P[facei],
                    pP_M[facei],
                    pR[facei],
                    pRT[facei],
                    rho_.value(),
                    pLmabda[facei],
                    pUp[facei],
                    pUs[facei]

                );
                //                 pFlux_ -> evaluateFlux
                // (
                //     pLmFlux[facei],
                //     pFFlux[facei],
                //     plmRight[facei],
                //     plmLeft[facei],
                //     pFRight[facei],
                //     pFLeft[facei],
                //     pPRight[facei],
                //     pPLeft[facei],
                //     pR[facei],
                //     pRT[facei],
                //     rho_.value(),
                //     pLmabda[facei],
                //     pUp[facei],
                //     pUs[facei]

                // );
            }

        }

    }


    // Update boundary field and values
    lm_b_.correctBoundaryConditions();
    t_b_.correctBoundaryConditions();

    forAll (lmFlux_.boundaryField(), patchi)
    {
        const fvPatch& curPatch = lm_.boundaryField()[patchi].patch();

        // Fluxes
#ifdef OPENFOAM_NOT_EXTEND
        fvsPatchVectorField& pLmFlux = lmFlux_.boundaryFieldRef()[patchi];
        fvsPatchTensorField& pFFlux = FFlux_.boundaryFieldRef()[patchi];
#else
        fvsPatchVectorField& pLmFlux = lmFlux_.boundaryField()[patchi];
        fvsPatchTensorField& pFFlux = FFlux_.boundaryField()[patchi];
#endif

    const labelList& fc =
                        lm_.boundaryField()[patchi].patch().faceCells();
                
        // const tensorField pPLeft  =  P_.boundaryField()[patchi].patchInternalField();
        const fvsPatchTensorField& pP_M = P_M_.boundaryField()[patchi];



        if(!curPatch.coupled())
        {
            forAll(curPatch, facei)
            {                
#ifdef OPENFOAM_NOT_EXTEND
                t_M_.boundaryFieldRef()[patchi][facei] =
                                      pP_M[facei] & N_.boundaryField()[patchi][facei];
                // t_M_.boundaryFieldRef()[patchi][facei] =
                //                       pPLeft[facei] & N_.boundaryField()[patchi][facei];
                                      
#else
                t_M_.boundaryField()[patchi][facei] =
                                      pP_M[facei] & N_.boundaryField()[patchi][facei];
                // t_M_.boundaryField()[patchi][facei] =
                //                       pPLeft[facei] & N_.boundaryField()[patchi][facei];

#endif


                pLmFlux[facei] =  
                    t_b_.boundaryField()[patchi][facei];

                pFFlux[facei] =  (1/rho_.value())*
                    ( lm_b_.boundaryField()[patchi][facei] * N_.boundaryField()[patchi][facei]);
        
            }
        }
    }
}


void Foam::solidNumericFlux::setRotationalMatrix()
{
    
// Loop through all faces in the mesh


#ifdef OPENFOAM_NOT_EXTEND
    const labelList& owner = mesh_.owner();
#else
    const auto& owner = mesh_.owner();
#endif


    forAll(owner, facei) 
    {
        // // Initialize an arbitrary vector
        vector arbitraryVector(1, 0, 0); // Default arbitrary vector along x-axis

        // Ensure the arbitrary vector is not parallel to the face normal N_[facei]
        if (mag(arbitraryVector & N_[facei]) > 0.999) // If nearly aligned
        {
            arbitraryVector = vector(0, 1, 0); // Switch to y-axis
            if (mag(arbitraryVector & N_[facei]) > 0.999) // If still nearly aligned
            {
                arbitraryVector = vector(0, 0, 1); // Use z-axis as a last resort
            }
        }

        // Step 3: Project arbitraryVector onto the plane tangent to N_[facei]
        T1_[facei] = arbitraryVector - (arbitraryVector & N_[facei]) * N_[facei];

        // Find an orthogonal vector to N_[facei]
        // T1_[facei] = op.findOrthogonal(N_[facei]);

        // Normalize T1_ to make it a unit vector
        scalar T1_Mag = mag(T1_[facei]);
        if (T1_Mag > SMALL) // Avoid division by zero
        {
            T1_[facei] /= T1_Mag;
        }
        else
        {
            FatalErrorInFunction << "Zero-length T1_ encountered at face " << facei << abort(FatalError);
        }

        // Second tangential vector (T2_): Cross product of N and T1_
        T2_[facei] = N_[facei] ^ T1_[facei];

        // Normalize T2_ to make it a unit vector
        scalar T2_Mag = mag(T2_[facei]);
        if (T2_Mag > SMALL) // Avoid division by zero
        {
            T2_[facei] /= T2_Mag;
        }
        else
        {
            FatalErrorInFunction << "Zero-length T2_ encountered at face " << facei << abort(FatalError);
        }

        // Verify orthogonality of T1_, T2_, and N
        scalar tolerance = 1e-6;

        if (mag(N_[facei] & T1_[facei]) > tolerance) 
        {
            FatalErrorInFunction << "N and T1_ are not orthogonal at face " << facei 
                                << ". Dot product: " << (N_[facei] & T1_[facei]) << abort(FatalError);
        }

        if (mag(N_[facei] & T2_[facei]) > tolerance) 
        {
            FatalErrorInFunction << "N and T2_ are not orthogonal at face " << facei 
                                << ". Dot product: " << (N_[facei] & T2_[facei]) << abort(FatalError);
        }

        if (mag(T1_[facei] & T2_[facei]) > tolerance) 
        {
            FatalErrorInFunction << "T1_ and T2_ are not orthogonal at face " << facei 
                                << ". Dot product: " << (T1_[facei] & T2_[facei]) << abort(FatalError);
        }

        // Now T1_ and T2_ are guaranteed orthogonal to N and each other


            //construct the rotational matrix
            R_[facei].xx() = N_[facei].x();
            R_[facei].xy() = N_[facei].y();
            R_[facei].xz() = N_[facei].z();
            R_[facei].yx() = T1_[facei].x();
            R_[facei].yy() = T1_[facei].y();
            R_[facei].yz() = T1_[facei].z();
            R_[facei].zx() = T2_[facei].x();
            R_[facei].zy() = T2_[facei].y();
            R_[facei].zz() = T2_[facei].z();

    }

// Coupled boundaries
forAll (lm_.boundaryField(), patchI)
{
    // if (mesh_.boundaryField()[patchI].coupled())
    // {
        // const fvPatch& curPatch = lm_.boundaryField()[patchI].patch();

        const labelList& fc =
            lm_.boundaryField()[patchI].patch().faceCells();
        
        // const fvPatch& curPatch = lm_.boundaryField()[patchi].patch();
       
       

#ifdef OPENFOAM_NOT_EXTEND
        const fvsPatchVectorField& pN = N_.boundaryField()[patchI];
         fvsPatchVectorField& pT1 = T1_.boundaryFieldRef()[patchI];
         fvsPatchVectorField& pT2 = T2_.boundaryFieldRef()[patchI];
         fvsPatchTensorField& pR = R_.boundaryFieldRef()[patchI];
#else
        const fvsPatchVectorField& pN = N_.boundaryField()[patchI];
         fvsPatchVectorField& pT1 = T1_.boundaryField()[patchI];
         fvsPatchVectorField& pT2 = T2_.boundaryField()[patchI];
         fvsPatchTensorField& pR = R_.boundaryField()[patchI];
#endif


        forAll (fc, facei)
        {
            // // Initialize an arbitrary vector
            vector arbitraryVector(1, 0, 0); // Default arbitrary vector along x-axis

            // Ensure the arbitrary vector is not parallel to the face normal N_[facei]
            if (mag(arbitraryVector & pN[facei]) > 0.999) // If nearly aligned
            {
                arbitraryVector = vector(0, 1, 0); // Switch to y-axis
                if (mag(arbitraryVector & pN[facei]) > 0.999) // If still nearly aligned
                {
                    arbitraryVector = vector(0, 0, 1); // Use z-axis as a last resort
                }
            }

            // Step 3: Project arbitraryVector onto the plane tangent to N_[facei]
            pT1[facei] = arbitraryVector - (arbitraryVector & pN[facei]) * pN[facei];


            // Normalize T1_ to make it a unit vector
            scalar T1_Mag = mag(pT1[facei]);
            if (T1_Mag > SMALL) // Avoid division by zero
            {
                pT1[facei] /= T1_Mag;
            }
            else
            {
                FatalErrorInFunction << "Zero-length T1_ encountered at boundary face " << facei << abort(FatalError);
            }

                // Info<< "face " << facei
                // << " |N| = " << mag(pN[facei]) << endl;

            // Second tangential vector (T2_): Cross product of N and T1_
            pT2[facei] = pN[facei] ^ pT1[facei];

            // Normalize T2_ to make it a unit vector
            scalar T2_Mag = mag(pT2[facei]);
            if (T2_Mag > SMALL) // Avoid division by zero
            {
                pT2[facei] /= T2_Mag;
            }
            else
            {
                FatalErrorInFunction << "Zero-length T2_ encountered at boundary face " << facei << abort(FatalError);
            }

            // Verify orthogonality of T1_, T2_, and N
            scalar tolerance = 1e-6;

            if (mag(pN[facei] & pT1[facei]) > tolerance) 
            {
                FatalErrorInFunction << "N and T1_ are not orthogonal at face " << facei 
                                    << ". Dot product: " << (pN[facei] & pT1[facei]) << abort(FatalError);
            }

            if (mag(pN[facei] & pT2[facei]) > tolerance) 
            {
                FatalErrorInFunction << "N and T2_ are not orthogonal at face " << facei 
                                    << ". Dot product: " << (pN[facei] & pT2[facei]) << abort(FatalError);
            }

            if (mag(pT1[facei] & pT2[facei]) > tolerance) 
            {
                FatalErrorInFunction << "T1_ and T2_ are not orthogonal at face " << facei 
                                    << ". Dot product: " << (pT1[facei] & pT2[facei]) << abort(FatalError);
            }

            // Now T1_ and T2_ are guaranteed orthogonal to N and each other


            //construct the rotational matrix
            pR[facei].xx() = pN[facei].x();
            pR[facei].xy() = pN[facei].y();
            pR[facei].xz() = pN[facei].z();
            pR[facei].yx() = pT1[facei].x();
            pR[facei].yy() = pT1[facei].y();
            pR[facei].yz() = pT1[facei].z();
            pR[facei].zx() = pT2[facei].x();
            pR[facei].zy() = pT2[facei].y();
            pR[facei].zz() = pT2[facei].z();
        }
    // }
}





RTranspos_ = R_.T();



}


// ************************************************************************* //
