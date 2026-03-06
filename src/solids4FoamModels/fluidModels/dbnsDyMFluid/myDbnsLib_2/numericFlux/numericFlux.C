/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     5.0
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
    q_
    (
        "q",
        dimensionSet(0, 2, -2, 0, 0),
        0
    ),
    pInf_
    (
        "p_inf",
        p_.dimensions(),
        0
    ),
    gamma_
    (
        "gamma",
        dimless,
        0
    ),
    cv_
    (
        "Cv",
        dimensionSet(0, 2, -2, -1, 0, 0, 0),
        0
    ),
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
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("fp", dimless, 1.0)
    ),

    fp1_
    {
        IOobject
        (
            "fp1",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("fp1", dimPressure/dimPressure, 1.0)
    }
{
   // Read thermodynamic properties from thermo dictionary
//    word firstToken;
//    thermo_.lookup("fluidProperties") >> firstToken >> q_.value() 
//    >> pInf_.value() >> gamma_.value() >> cv_.value();
 word thermoTypeName;
            
                IOdictionary thermoDict
                (
                    IOobject
                    (
                        "thermophysicalProperties",
                        mesh_.time().constant(),
                        mesh_,
                        IOobject::MUST_READ_IF_MODIFIED,
                        IOobject::NO_WRITE
                    )
                );

                thermoDict.lookup("thermoType") >> thermoTypeName;
            
            if (thermoTypeName == "externalStiffenedGasThermo")
            {
                const dictionary& d = thermoDict.subDict("stiffenedGasCoeffs");
                d.lookup("pInf") >> pInf_.value();
                d.lookup("gamma") >> gamma_.value();
                d.lookup("Cv") >> cv_.value();
                d.lookup("q") >> q_.value();
            }
            else
            {
                FatalErrorIn("hllcSGLMFlux")
                << "thermo Type Name = " << thermoTypeName << nl
                << " it should be externalStiffenedGasThermo " << nl
                << exit(FatalError);

                //   pinf =0;
                //  gamma = 1.4;
                //  Cv = 718;
                //  q =0;
            }
            Info<< "Fluid Properties:"<< nl
                << "pInf= "<< pInf_.value() << nl
                << "gamma= "<< gamma_.value() << nl
                << "Cv= "<< cv_.value() << nl
                << "q= "<< q_.value() << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Flux, class Limiter>
void Foam::numericFlux<Flux, Limiter>::computeFlux()
{
    // Get face-to-cell addressing: face area point from owner to neighbour
    const unallocLabelList& owner = this->mesh().owner();
    const unallocLabelList& neighbour = this->mesh().neighbour();

    // Get the face area vector
    const surfaceVectorField& Sf = this->mesh().Sf();
    const surfaceScalarField& magSf = this->mesh().magSf();

    const volVectorField& cellCentre = this->mesh().C();
    const surfaceVectorField& faceCentre = this->mesh().Cf();

    // ALE mesh velocity 
    surfaceScalarField mshPhi( meshPhi() ); 

    // Get gradients
    // Coupled patch update on gradients moved into gradScheme.C
    // HJ, 22/Apr/2016;

    // Changed return type for gradient cacheing.  HJ, 22/Apr/2016
    const tmp<volVectorField> tgradP = fvc::grad(p_);
    const volVectorField& gradP = tgradP();

    const tmp<volTensorField> tgradU = fvc::grad(U_);
    const volTensorField& gradU = tgradU();

    const tmp<volVectorField> tgradT = fvc::grad(T_);
    const volVectorField& gradT = tgradT();

    MDLimiter<scalar, Limiter> scalarPLimiter
    (
        this->p_,
        gradP
    );

    MDLimiter<vector, Limiter> vectorULimiter
    (
        this->U_,
        gradU
    );

    MDLimiter<scalar, Limiter> scalarTLimiter
    (
        this->T_,
        gradT
    );

    // Get limiters
    const volScalarField& pLimiter = scalarPLimiter.phiLimiter();
    const volVectorField& ULimiter = vectorULimiter.phiLimiter();
    const volScalarField& TLimiter = scalarTLimiter.phiLimiter();

    //!---------------------------------------------
    // evalautePressureSensor();

    //!---------------------------------------------


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
            p_[own] + pLimiter[own]*(deltaRLeft & gradP[own]),
            p_[nei] + pLimiter[nei]*(deltaRRight & gradP[nei]),
            U_[own] + cmptMultiply(ULimiter[own], (deltaRLeft & gradU[own])),
            U_[nei] + cmptMultiply(ULimiter[nei], (deltaRRight & gradU[nei])),
            T_[own] + TLimiter[own]*(deltaRLeft & gradT[own]),
            T_[nei] + TLimiter[nei]*(deltaRRight & gradT[nei]),
            q_.value(),
            pInf_.value(),
            gamma_.value(),
            cv_.value(),
            Sf[faceI],
            magSf[faceI],
            mshPhi[faceI],
            fp1_[own],
            fp1_[nei]
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

        // const scalarField& pcv = cv.boundaryField()[patchi];
        // const scalarField& pR = R.boundaryField()[patchi];

        // Gradients
        const fvPatchVectorField& pGradP = gradP.boundaryField()[patchi];
        const fvPatchTensorField& pGradU = gradU.boundaryField()[patchi];
        const fvPatchVectorField& pGradT = gradT.boundaryField()[patchi];

        // Limiters
        const fvPatchScalarField& pPatchLim = pLimiter.boundaryField()[patchi];
        const fvPatchVectorField& UPatchLim = ULimiter.boundaryField()[patchi];
        const fvPatchScalarField& TPatchLim = TLimiter.boundaryField()[patchi];

        // Face areas
        const fvsPatchVectorField& pSf = Sf.boundaryField()[patchi];
        const fvsPatchScalarField& pMagSf = magSf.boundaryField()[patchi];
        const fvsPatchScalarField& pMshPhi = mshPhi.boundaryField()[patchi];


        const scalarField& pfp1  = fp1_.boundaryField()[patchi];


        if (pp.coupled())
        {
            // Coupled patch
            const scalarField ppLeft  =
                p_.boundaryField()[patchi].patchInternalField();

            const scalarField ppRight =
                p_.boundaryField()[patchi].patchNeighbourField();

            const vectorField pULeft  =
                U_.boundaryField()[patchi].patchInternalField();

            const vectorField pURight =
                U_.boundaryField()[patchi].patchNeighbourField();

            const scalarField pTLeft  =
                T_.boundaryField()[patchi].patchInternalField();

            const scalarField pTRight =
                T_.boundaryField()[patchi].patchNeighbourField();

            // Gradients
            const vectorField pgradPLeft = pGradP.patchInternalField();
            const vectorField pgradPRight = pGradP.patchNeighbourField();

            const tensorField pgradULeft = pGradU.patchInternalField();
            const tensorField pgradURight = pGradU.patchNeighbourField();

            const vectorField pgradTLeft = pGradT.patchInternalField();
            const vectorField pgradTRight = pGradT.patchNeighbourField();

            // Geometry: call the raw cell-to-face vector by calling
            // the base patch (cell-to-face) delta coefficient
            // Work out the right delta from the cell-to-cell delta
            // across the coupled patch and left delta
            vectorField pDeltaRLeft = curPatch.fvPatch::delta();
            vectorField pDdeltaRRight = pDeltaRLeft - curPatch.delta();

            // Limiters

            const scalarField ppLimiterLeft = pPatchLim.patchInternalField();
            const scalarField ppLimiterRight = pPatchLim.patchNeighbourField();

            const vectorField pULimiterLeft = UPatchLim.patchInternalField();
            const vectorField pULimiterRight = UPatchLim.patchNeighbourField();

            const scalarField pTLimiterLeft = TPatchLim.patchInternalField();
            const scalarField pTLimiterRight = TPatchLim.patchNeighbourField();

            forAll (pp, facei)
            {
                Flux::evaluateFlux
                (
                    pRhoFlux[facei],
                    pRhoUFlux[facei],
                    pRhoEFlux[facei],

                    ppLeft[facei]
                + ppLimiterLeft[facei]*
                    (pDeltaRLeft[facei] & pgradPLeft[facei]),

                    ppRight[facei]
                + ppLimiterRight[facei]*
                    (pDdeltaRRight[facei] & pgradPRight[facei]),

                    pULeft[facei]
                + cmptMultiply
                    (
                        pULimiterLeft[facei],
                        pDeltaRLeft[facei] & pgradULeft[facei]
                    ),

                    pURight[facei]
                + cmptMultiply
                    (
                        pULimiterRight[facei],
                        pDdeltaRRight[facei] & pgradURight[facei]
                    ),

                    pTLeft[facei]
                + pTLimiterLeft[facei]*
                    (pDeltaRLeft[facei] & pgradTLeft[facei]),

                    pTRight[facei]
                + pTLimiterRight[facei]*
                    (pDdeltaRRight[facei] & pgradTRight[facei]),

                    q_.value(),
                    pInf_.value(),
                    gamma_.value(),
                    cv_.value(),
                    pSf[facei],
                    pMagSf[facei],
                    pMshPhi[facei],
                    pfp1[facei],
                    pfp1[facei]
                    
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
                    q_.value(),
                    pInf_.value(),
                    gamma_.value(),
                    cv_.value(),
                    pSf[facei],
                    pMagSf[facei],
                    pMshPhi[facei],
                    pfp1[facei],
                    pfp1[facei]
                );
            }
        } 
    }
}


// ************************************************************************* //
template<class Flux, class Limiter>
void Foam::numericFlux<Flux, Limiter>::evalautePressureSensor
(
    // const surfaceScalarField& p_neg,
    // const surfaceScalarField& p_pos
)
{
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();

    // ----------------------------------------------------
    // Step 0: Initialize cell sensor to large value
    // ----------------------------------------------------
    fp1_ = scalar(GREAT);

    // ----------------------------------------------------
    // Step 1: Compute face sensor (internal faces)
    // ----------------------------------------------------
    forAll(owner, faceI)
    {

        const label own = owner[faceI];
        const label nei = neighbour[faceI];

        scalar ratio =
            min
            (
                p_[own] / (p_[nei] + VSMALL),
                p_[nei] / (p_[own] + VSMALL)
            );

        fp_[faceI] = pow3(ratio);
    }

    // ----------------------------------------------------
    // Step 2: Compute face sensor (boundary faces)
    // ----------------------------------------------------
    forAll(mesh_.boundary(), patchI)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];

        if (mesh_.boundary()[patchI].type() == "empty")
        {
            continue;
        }
        
        const fvPatchScalarField& pp = p_.boundaryField()[patchI];

        forAll(patch, facei)
        {

            scalar ratio = 1;
                // min
                // (
                //     p_neg.boundaryField()[patchI][facei]
                //   / (p_pos.boundaryField()[patchI][facei] + VSMALL),

                //     p_pos.boundaryField()[patchI][facei]
                //   / (p_neg.boundaryField()[patchI][facei] + VSMALL)
                // );

            fp_.boundaryField()[patchI][facei] = pow3(ratio);
        }
    }

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

