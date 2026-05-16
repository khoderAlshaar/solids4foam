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
#include "numericFlux.H"
#include "MDLimiter.H"
#include "IOstreams.H"
#include "tmp.H" 
// #include "directionInterpolate.H"
#include "slipFvPatchFields.H"
#include "movingWallVelocityFvPatchVectorField.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
template<class Flux, class Limiter>
Foam::numericFlux<Flux, Limiter>::numericFlux
(
    const volScalarField& p,
    const volVectorField& U,
    const volScalarField& T,
    basicThermo& thermo
)
:
    numericFluxBase<Flux>(p.mesh()),
    mesh_(p.mesh()),
    p_(p),
    U_(U),
    T_(T),
    thermo_(thermo),
    rhoFlux_
    (
        IOobject
        (
            "phi",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (linearInterpolate(thermo_.rho()*U_) & this->mesh().Sf())
    ),
    rhoUFlux_
    (
        IOobject
        (
            "rhoUFlux",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rhoFlux_*linearInterpolate(U_)
    ),
    rhoEFlux_
    (
        IOobject
        (
            "rhoEFlux",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rhoFlux_*linearInterpolate(thermo.Cv()*T_ + 0.5*magSqr(U_))
    ),
    fp_
    (
        IOobject
        (
            "fp",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("fp", dimless, 1.0)
    ),
    fp1_
    (
        IOobject
        (
            "fp1",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("fp1", dimPressure/dimPressure, 1.0)
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
    )
{
       setRotationalMatrix();

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Flux, class Limiter>
void Foam::numericFlux<Flux, Limiter>::computeFlux()
{
    // Get face-to-cell addressing: face area point from owner to neighbour
    const auto& owner = this->mesh().owner();
    const auto& neighbour = this->mesh().neighbour();

    // Get the face area vector
    const surfaceVectorField& Sf = this->mesh().Sf();
    const surfaceScalarField& magSf = this->mesh().magSf();

    const volVectorField& cellCentre = this->mesh().C();
    const surfaceVectorField& faceCentre = this->mesh().Cf();

    // ALE mesh velocity 
    surfaceScalarField mshPhi( meshPhi() ); 

    // Thermodynamics
    const volScalarField Cv = thermo_.Cv();
    const volScalarField R  = thermo_.Cp() - Cv;

    // Get gradients
    // Coupled patch update on gradients moved into gradScheme.C
    // HJ, 22/Apr/2016;

    // Changed return type for gradient cacheing.  HJ, 22/Apr/2016
    // const tmp<volVectorField> tgradP = fvc::grad(p_);
    // const volVectorField& gradP = tgradP();

    // const tmp<volTensorField> tgradU = fvc::grad(U_);
    // const volTensorField& gradU = tgradU();

    // const tmp<volVectorField> tgradT = fvc::grad(T_);
    // const volVectorField& gradT = tgradT();

    // MDLimiter<scalar, Limiter> scalarPLimiter
    // (
    //     this->p_,
    //     gradP
    // );

    // MDLimiter<vector, Limiter> vectorULimiter
    // (
    //     this->U_,
    //     gradU
    // );

    // MDLimiter<scalar, Limiter> scalarTLimiter
    // (
    //     this->T_,
    //     gradT
    // );

    // // Get limiters
    // const volScalarField& pLimiter = scalarPLimiter.phiLimiter();
    // const volVectorField& ULimiter = vectorULimiter.phiLimiter();
    // const volScalarField& TLimiter = scalarTLimiter.phiLimiter();

    surfaceScalarField pos_(IOobject("pos", mesh_), mesh_, dimensionedScalar("one", dimless, 1.0));
    surfaceScalarField neg_(IOobject("neg", mesh_), mesh_, dimensionedScalar("minusOne", dimless, -1.0));

    surfaceScalarField p_pos( fvc::interpolate(p_, pos_,"reconstruct(p)") );
    surfaceScalarField p_neg( fvc::interpolate(p_, neg_,"reconstruct(p)") );

    surfaceVectorField U_pos( fvc::interpolate(U_, pos_,"reconstruct(U)") );
    surfaceVectorField U_neg( fvc::interpolate(U_, neg_,"reconstruct(U)") );

    surfaceScalarField T_pos( fvc::interpolate(T_, pos_,"reconstruct(T)") );
    surfaceScalarField T_neg( fvc::interpolate(T_, neg_,"reconstruct(T)") );


    evalautePressureSensor(p_);


    // Calculate fluxes at internal faces
    forAll (owner, faceI)
    {
        const label own = owner[faceI];
        const label nei = neighbour[faceI];

        const vector deltaRLeft = faceCentre[faceI] - cellCentre[own];
        const vector deltaRRight = faceCentre[faceI] - cellCentre[nei];

        // calculate fluxes with reconstructed primitive variables at faces
        Flux::evaluateFlux
        (
            rhoFlux_[faceI],
            rhoUFlux_[faceI],
            rhoEFlux_[faceI],
            // p_[own]  + pLimiter[own]*(deltaRLeft & gradP[own]),
            // p_[nei]  + pLimiter[nei]*(deltaRRight & gradP[nei]),
            // U_[own]  + cmptMultiply(ULimiter[own], (deltaRLeft & gradU[own])),
            // U_[nei]  + cmptMultiply(ULimiter[nei], (deltaRRight & gradU[nei])),
            // T_[own]  + TLimiter[own]*(deltaRLeft & gradT[own]),
            // T_[nei]  + TLimiter[nei]*(deltaRRight & gradT[nei]),
            p_pos[faceI],  p_neg[faceI],
            U_pos[faceI],  U_neg[faceI],
            T_pos[faceI],  T_neg[faceI],
            R[own],
            R[nei],
            Cv[own],
            Cv[nei],
            Sf[faceI],
            magSf[faceI],
            mshPhi[faceI],
            fp1_[own],
            fp1_[nei],
            R_[faceI],
            RTranspos_[faceI]
        );
    }

    // Update boundary field and values
    forAll (rhoFlux_.boundaryField(), patchi)
    {
        const fvPatch& curPatch = p_.boundaryField()[patchi].patch();

        // Fluxes
        fvsPatchScalarField& pRhoFlux  = rhoFlux_.boundaryField()[patchi];
        fvsPatchVectorField& pRhoUFlux = rhoUFlux_.boundaryField()[patchi];
        fvsPatchScalarField& pRhoEFlux = rhoEFlux_.boundaryField()[patchi];

        // Patch fields
        const fvPatchScalarField& pp = p_.boundaryField()[patchi];
        const vectorField& pU = U_.boundaryField()[patchi];
        const scalarField& pT = T_.boundaryField()[patchi];

        const scalarField& pCv = Cv.boundaryField()[patchi];
        const scalarField& pR = R.boundaryField()[patchi];

        const scalarField& pfp1  = fp1_.boundaryField()[patchi];


        // Gradients
        // const fvPatchVectorField& pGradP = gradP.boundaryField()[patchi];
        // const fvPatchTensorField& pGradU = gradU.boundaryField()[patchi];
        // const fvPatchVectorField& pGradT = gradT.boundaryField()[patchi];

        // // Limiters
        // const fvPatchScalarField& pPatchLim = pLimiter.boundaryField()[patchi];
        // const fvPatchVectorField& UPatchLim = ULimiter.boundaryField()[patchi];
        // const fvPatchScalarField& TPatchLim = TLimiter.boundaryField()[patchi];

        // Face areas
        const fvsPatchVectorField& pSf = Sf.boundaryField()[patchi];
        const fvsPatchScalarField& pMagSf = magSf.boundaryField()[patchi];
        const fvsPatchScalarField& pMshPhi = mshPhi.boundaryField()[patchi];

            // Face areas
        const fvsPatchTensorField& pRR = R_.boundaryField()[patchi];
        const fvsPatchTensorField& pRT = RTranspos_.boundaryField()[patchi];

        const scalarField ppLeft  =  p_.boundaryField()[patchi].patchInternalField();
        const scalarField pTLeft  =  T_.boundaryField()[patchi].patchInternalField();
        const vectorField pULeft  =  U_.boundaryField()[patchi].patchInternalField();
        
        if (curPatch.coupled())
        {
            //             // Coupled patch
            // const scalarField ppLeft  =
            //     p_.boundaryField()[patchi].patchInternalField();

            // const scalarField ppRight =
            //     p_.boundaryField()[patchi].patchNeighbourField();

            // const vectorField pULeft  =
            //     U_.boundaryField()[patchi].patchInternalField();

            // const vectorField pURight =
            //     U_.boundaryField()[patchi].patchNeighbourField();

            // const scalarField pTLeft  =
            //     T_.boundaryField()[patchi].patchInternalField();

            // const scalarField pTRight =
            //     T_.boundaryField()[patchi].patchNeighbourField();
            // // Gradients
            // const vectorField pgradPLeft = pGradP.patchInternalField();
            // const vectorField pgradPRight = pGradP.patchNeighbourField();

            // const tensorField pgradULeft = pGradU.patchInternalField();
            // const tensorField pgradURight = pGradU.patchNeighbourField();

            // const vectorField pgradTLeft = pGradT.patchInternalField();
            // const vectorField pgradTRight = pGradT.patchNeighbourField();

            // Geometry: call the raw cell-to-face vector by calling
            // the base patch (cell-to-face) delta coefficient
            // Work out the right delta from the cell-to-cell delta
            // across the coupled patch and left delta
            vectorField pDeltaRLeft = curPatch.fvPatch::delta();
            vectorField pDdeltaRRight = pDeltaRLeft - curPatch.delta();

            // Limiters

            // const scalarField ppLimiterLeft = pPatchLim.patchInternalField();
            // const scalarField ppLimiterRight = pPatchLim.patchNeighbourField();

            // const vectorField pULimiterLeft = UPatchLim.patchInternalField();
            // const vectorField pULimiterRight = UPatchLim.patchNeighbourField();

            // const scalarField pTLimiterLeft = TPatchLim.patchInternalField();
            // const scalarField pTLimiterRight = TPatchLim.patchNeighbourField();

            // Patch fields
            const fvsPatchScalarField& pp_pos = p_pos.boundaryField()[patchi];
            const fvsPatchVectorField& pU_pos = U_pos.boundaryField()[patchi];
            const fvsPatchScalarField& pT_pos = T_pos.boundaryField()[patchi];
            
            const fvsPatchScalarField& pp_neg = p_neg.boundaryField()[patchi];
            const fvsPatchVectorField& pU_neg = U_neg.boundaryField()[patchi];
            const fvsPatchScalarField& pT_neg = T_neg.boundaryField()[patchi];
            
            forAll (curPatch, facei)
            {
                Flux::evaluateFlux
                (
                    pRhoFlux[facei],
                    pRhoUFlux[facei],
                    pRhoEFlux[facei],

                    pp_pos[facei],  pp_neg[facei],
                    pU_pos[facei],  pU_neg[facei],
                    pT_pos[facei],  pT_neg[facei],
                //                         ppLeft[facei]
                // + ppLimiterLeft[facei]*
                //     (pDeltaRLeft[facei] & pgradPLeft[facei]),

                //     ppRight[facei]
                // + ppLimiterRight[facei]*
                //     (pDdeltaRRight[facei] & pgradPRight[facei]),

                //     pULeft[facei]
                // + cmptMultiply
                //     (
                //         pULimiterLeft[facei],
                //         pDeltaRLeft[facei] & pgradULeft[facei]
                //     ),

                //     pURight[facei]
                // + cmptMultiply
                //     (
                //         pULimiterRight[facei],
                //         pDdeltaRRight[facei] & pgradURight[facei]
                //     ),

                //     pTLeft[facei]
                // + pTLimiterLeft[facei]*
                //     (pDeltaRLeft[facei] & pgradTLeft[facei]),

                //     pTRight[facei]
                // + pTLimiterRight[facei]*
                //     (pDdeltaRRight[facei] & pgradTRight[facei]),


                    pR[facei],  pR[facei],
                    pCv[facei], pCv[facei],
                    pSf[facei],
                    pMagSf[facei],
                    pMshPhi[facei],
                    pfp1[facei],
                    pfp1[facei],
                    pRR[facei],
                    pRT[facei]
                );
            }
        }
        else if (curPatch.type() == "wall")
        // if (curPatch.type() == "wall")
        {
                        // Gradients
            // const vectorField pgradPLeft = pGradP.patchInternalField();

            // const tensorField pgradULeft = pGradU.patchInternalField();

            // const vectorField pgradTLeft = pGradT.patchInternalField();

            vectorField pDeltaRLeft = curPatch.fvPatch::delta();

                        // Limiters

            // const scalarField ppLimiterLeft = pPatchLim.patchInternalField();

            // const vectorField pULimiterLeft = UPatchLim.patchInternalField();

            // const scalarField pTLimiterLeft = TPatchLim.patchInternalField();

            // Info << "The BCS is Wall"<< endl;
            const fvPatchScalarField& pp = p_.boundaryField()[patchi];
            const fvPatchVectorField& pU = U_.boundaryField()[patchi];
            const scalarField& pT = T_.boundaryField()[patchi];

            // Patch fields
            const fvsPatchScalarField& pp_pos = p_pos.boundaryField()[patchi];
            const fvsPatchVectorField& pU_pos = U_pos.boundaryField()[patchi];
            const fvsPatchScalarField& pT_pos = T_pos.boundaryField()[patchi];

            forAll (pp, facei)
            {
                const label own = owner[facei];
                const label nei = neighbour[facei];
                const vector n = pSf[facei] / pMagSf[facei];

                // Reconstructed left state
                const vector ULeft =   pU_pos[facei];
                // const vector ULeft =
                //     pULeft[facei]
                // + cmptMultiply(pULimiterLeft[facei],
                //                 pDeltaRLeft[facei] & pgradULeft[facei]);

                // Wall velocity from BC
                const vector Uwall = pU[facei];

                vector URight;

                if (isA<slipFvPatchVectorField>(pU) ||isA<movingWallVelocityFvPatchVectorField>(pU))
                {
                    // -------------------------
                    // SLIP WALL (moving or not)
                    // -------------------------
                    const vector Urel = ULeft - Uwall;
                    URight = Uwall + (Urel - 2.0*(Urel & n)*n);
                }
                else
                {
                    // -------------------------
                    // NO-SLIP / MOVING WALL
                    // -------------------------
                    URight = 2.0*Uwall - ULeft;
                }

                // Calculate fluxes
                Flux::evaluateFlux
                (
                    pRhoFlux[facei],
                    pRhoUFlux[facei],
                    pRhoEFlux[facei],
                    // ppLeft[facei] + ppLimiterLeft[facei]*(pDeltaRLeft[facei] & pgradPLeft[facei]),
                    // ppLeft[facei] + ppLimiterLeft[facei]*(pDeltaRLeft[facei] & pgradPLeft[facei]),
                    // ULeft,
                    // URight,
                    // pTLeft[facei] + pTLimiterLeft[facei]*(pDeltaRLeft[facei] & pgradTLeft[facei]),
                    // pTLeft[facei] + pTLimiterLeft[facei]*(pDeltaRLeft[facei] & pgradTLeft[facei]),

                    pp_pos[facei],
                    pp_pos[facei],
                    ULeft,
                    URight,
                    pT_pos[facei],
                    pT_pos[facei],
                    // ppLeft[facei],
                    // pp[facei],
                    // pULeft[facei],
                    // pU[facei],
                    // pTLeft[facei],
                    // pT[facei],
                    pR[facei],
                    pR[facei],
                    pCv[facei],
                    pCv[facei],
                    pSf[facei],
                    pMagSf[facei],
                    pMshPhi[facei],
                    pfp1[facei],
                    pfp1[facei],
                    pRR[facei],
                    pRT[facei]
                );
            }
        }
        else
        {
            forAll (pp, facei)
            {
                // Calculate fluxes
                Flux::evaluateFlux
                (
                    pRhoFlux[facei],
                    pRhoUFlux[facei],
                    pRhoEFlux[facei],
                    pp[facei],
                    pp[facei],
                    pU[facei],
                    pU[facei],
                    pT[facei],
                    pT[facei],
                    // ppLeft[facei],
                    // pp[facei],
                    // pULeft[facei],
                    // pU[facei],
                    // pTLeft[facei],
                    // pT[facei],
                    pR[facei],
                    pR[facei],
                    pCv[facei],
                    pCv[facei],
                    pSf[facei],
                    pMagSf[facei],
                    pMshPhi[facei],
                    pfp1[facei],
                    pfp1[facei],
                    pRR[facei],
                    pRT[facei]
                );
            }
        }
    }
}


// ************************************************************************* //
template<class Flux, class Limiter>
Foam::tmp<Foam::surfaceScalarField> numericFlux<Flux, Limiter>::meshPhi() const
{
    if (mesh_.moving()) 
    {
        return  fvc::meshPhi(U_);
    } 

    return tmp<surfaceScalarField>
        (
            new surfaceScalarField
            (
                IOobject
                (
                "meshPhi",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
                ),
                mesh_,
                dimensionedScalar("0", dimVolume/dimTime, 0.0)
            )
        );
    
}



// ************************************************************************* //


template<class Flux, class Limiter>
void Foam::numericFlux<Flux, Limiter>::evalautePressureSensor
(
    const volScalarField& p
)
{

    // surfaceScalarField p_neg
    // (
    //     IOobject
    //     (
    //         "p_neg",
    //         this->mesh().time().timeName(),
    //         this->mesh(),
    //         IOobject::NO_READ,
    //         IOobject::NO_WRITE
    //     ),
    //     p
    // );
    // surfaceScalarField p_pos
    // (
    //     IOobject
    //     (
    //         "p_pos",
    //         this->mesh().time().timeName(),
    //         this->mesh(),
    //         IOobject::NO_READ,
    //         IOobject::NO_WRITE
    //     ),
    //     p
    // );
    //     // Calculate fluxes at internal faces
    // forAll (owner, faceI)
    // {
    //     const label own = owner[faceI];
    //     const label nei = neighbour[faceI];

    //     // const vector deltaRLeft = faceCentre[faceI] - cellCentre[own];
    //     // const vector deltaRRight = faceCentre[faceI] - cellCentre[nei];

    //    p_neg[faceI] = p_[own] ;//+ pLimiter[own]*(deltaRLeft & gradP[own]),
    //    p_pos[faceI] = p_[nei] ;//+ pLimiter[nei]*(deltaRRight & gradP[nei]),
    // }
    
    const labelUList& owner = mesh().owner();
    const labelUList& neighbour = mesh().neighbour();

    // ----------------------------------------------------
    // Step 0: Initialize cell sensor to large value
    // ----------------------------------------------------
    fp1_ = scalar(1.0);

    // ----------------------------------------------------
    // Step 1: Compute face sensor (internal faces)
    // ----------------------------------------------------
    forAll(owner, faceI)
    {
        const label own = owner[faceI];
        const label nei = neighbour[faceI];
        scalar  p_neg = p_[own] ;
        scalar  p_pos = p_[nei] ;
        scalar ratio =
            min
            (
                p_neg / (p_pos + VSMALL),
                p_pos / (p_neg + VSMALL)
            );

        fp_[faceI] = pow3(ratio);
    }

    // ----------------------------------------------------
    // Step 2: Compute face sensor (boundary faces)
    // ----------------------------------------------------
    
    // forAll(mesh_.boundary(), patchI)
    // {
    //     const fvPatch& patch = mesh_.boundary()[patchI];

    //     const fvPatchScalarField& pp = p_.boundaryField()[patchI];


    //     if (mesh_.boundary()[patchI].type() == "empty")
    //     {
    //         continue;
    //     }
        
    //     forAll(patch, facei)
    //     {
    //                 scalar  p_neg[faceI] = p_[own] ;//+ pLimiter[own]*(deltaRLeft & gradP[own]),
    //     scalar  p_pos[faceI] = p_[nei] ;

    //         scalar ratio =
    //             min
    //             (
    //                 p_neg.boundaryField()[patchI][facei]
    //               / (p_pos.boundaryField()[patchI][facei] + VSMALL),

    //                 p_pos.boundaryField()[patchI][facei]
    //               / (p_neg.boundaryField()[patchI][facei] + VSMALL)
    //             );

    //         fp_.boundaryField()[patchI][facei] = pow3(ratio);
    //     }
    // }

    // ----------------------------------------------------
    // Step 3: Reduce face sensor → cell sensor (internal)
    // ----------------------------------------------------
    forAll(owner, faceI)
    {
        const label own = owner[faceI];
        const label nei = neighbour[faceI];

        const scalar val = fp_[faceI];

        fp1_[own] = min(fp1_[own], val);
        fp1_[nei] = min(fp1_[nei], val);
    }

    // ----------------------------------------------------
    // Step 4: Reduce face sensor → cell sensor (boundary)
    // ----------------------------------------------------
    forAll(mesh_.boundary(), patchI)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];

        if (mesh_.boundary()[patchI].type() == "empty")
        {
            continue;
        }

        const scalarField& fPatch = fp_.boundaryField()[patchI];
        const labelUList& faceCells = patch.faceCells();

        forAll(patch, facei)
        {
            label cellI = faceCells[facei];
            fp1_[cellI] = min(fp1_[cellI], fPatch[facei]);
        }


        scalarField& bField = fp1_.boundaryField()[patchI];
        // const labelUList& faceCells = patch.faceCells();

        forAll(patch, facei)
        {
            bField[facei] = fp1_[faceCells[facei]];
        }
    }

    // fp1_.correctBoundaryConditions();
}



template<class Flux, class Limiter>
void Foam::numericFlux<Flux, Limiter>::setRotationalMatrix()
{

// #ifdef OPENFOAM_NOT_EXTEND
//     const labelList& owner = mesh().owner();
// #else
    const auto& owner = mesh().owner();
// #endif

    const scalar tol = 1e-8;
    // Get the face area vector
    const surfaceVectorField& Sf = mesh().Sf();
    const surfaceScalarField& magSf = mesh().magSf();

    // =========================
    // Internal faces
    // =========================
    forAll(owner, facei)
    {
        vector n = Sf[facei]/magSf[facei];

        scalar magN = mag(n);

        if (magN < SMALL)
        {
            FatalErrorInFunction
                << "Zero normal at face " << facei
                << abort(FatalError);
        }

        // Normalize
        n /= magN;

        scalar n1 = n.x();
        scalar n2 = n.y();
        scalar n3 = n.z();

        // s = sign(n1), fallback = 1
        scalar s = (mag(n1) > SMALL) ? sign(n1) : 1.0;

        scalar denom = n1 + s;

        if (mag(denom) < SMALL)
        {
            FatalErrorInFunction
                << "Singular rotation matrix at face " << facei
                << " (n1 + s ≈ 0)"
                << abort(FatalError);
        }

        // =========================
        // Build rotation matrix
        // =========================
        tensor& Rf = R_[facei];

        Rf.xx() = n1;
        Rf.xy() = n2;
        Rf.xz() = n3;

        Rf.yx() = -n2;
        Rf.yy() = n1 + (n3*n3)/denom;
        Rf.yz() = -(n2*n3)/denom;

        Rf.zx() = -n3;
        Rf.zy() = -(n3*n2)/denom;
        Rf.zz() = n1 + (n2*n2)/denom;

        // =========================
        // DEBUG: Check R * R^T = I
        // =========================
        tensor Icheck = Rf & Rf.T();

        if
        (
            mag(Icheck.xx() - 1) > tol ||
            mag(Icheck.yy() - 1) > tol ||
            mag(Icheck.zz() - 1) > tol ||
            mag(Icheck.xy()) > tol ||
            mag(Icheck.xz()) > tol ||
            mag(Icheck.yx()) > tol ||
            mag(Icheck.yz()) > tol ||
            mag(Icheck.zx()) > tol ||
            mag(Icheck.zy()) > tol
        )
        {
            FatalErrorInFunction
                << "Rotation matrix is not orthogonal at face " << facei << nl
                << "R = " << Rf << nl
                << "R*R^T = " << Icheck << nl
                << abort(FatalError);
        }
    }

    // =========================
    // Boundary faces
    // =========================
    forAll(p_.boundaryField(), patchI)
    {
        const labelList& fc =
            p_.boundaryField()[patchI].patch().faceCells();

// #ifdef OPENFOAM_NOT_EXTEND
//         const fvsPatchVectorField& pN = N_.boundaryField()[patchI];
//         fvsPatchTensorField& pR = R_.boundaryFieldRef()[patchI];
// #else
        // Face areas
        const fvsPatchVectorField& pSf = Sf.boundaryField()[patchI];
        const fvsPatchScalarField& pMagSf = magSf.boundaryField()[patchI];
        fvsPatchTensorField& pR = R_.boundaryField()[patchI];
// #endif

        forAll(fc, facei)
        {
            vector n =  pSf[facei]/pMagSf[facei];
            scalar magN = mag(n);

            if (magN < SMALL)
            {
                FatalErrorInFunction
                    << "Zero normal at boundary face " << facei
                    << abort(FatalError);
            }

            n /= magN;

            scalar n1 = n.x();
            scalar n2 = n.y();
            scalar n3 = n.z();

            scalar s = (mag(n1) > SMALL) ? sign(n1) : 1.0;
            scalar denom = n1 + s;

            if (mag(denom) < SMALL)
            {
                FatalErrorInFunction
                    << "Singular rotation matrix at boundary face "
                    << facei
                    << abort(FatalError);
            }

            tensor& Rf = pR[facei];

            Rf.xx() = n1;
            Rf.xy() = n2;
            Rf.xz() = n3;

            Rf.yx() = -n2;
            Rf.yy() = n1 + (n3*n3)/denom;
            Rf.yz() = -(n2*n3)/denom;

            Rf.zx() = -n3;
            Rf.zy() = -(n3*n2)/denom;
            Rf.zz() = n1 + (n2*n2)/denom;

            // =========================
            // DEBUG CHECK
            // =========================
            tensor Icheck = Rf & Rf.T();

            if
            (
                mag(Icheck.xx() - 1) > tol ||
                mag(Icheck.yy() - 1) > tol ||
                mag(Icheck.zz() - 1) > tol ||
                mag(Icheck.xy()) > tol ||
                mag(Icheck.xz()) > tol ||
                mag(Icheck.yx()) > tol ||
                mag(Icheck.yz()) > tol ||
                mag(Icheck.zx()) > tol ||
                mag(Icheck.zy()) > tol
            )
            {
                FatalErrorInFunction
                    << "Rotation matrix NOT orthogonal at boundary face "
                    << facei << nl
                    << "R = " << Rf << nl
                    << "R*R^T = " << Icheck << nl
                    << abort(FatalError);
            }
        }
    }

    // =========================
    // Store transpose
    // =========================
    RTranspos_ = R_.T();
}