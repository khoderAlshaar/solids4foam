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
// #include "directionInterpolate.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::numericFlux::numericFlux
(
    const volScalarField& p,
    const volVectorField& U,
    const volScalarField& T,
    basicThermo& thermo
    // const MRFZoneList& MRF
)
:
    pFlux_( Foam::dbnsFlux::New( 
        p.mesh(), 
        p.mesh().thisDb().lookupObject<IOdictionary>("fvSchemes")) 
    ),
    mesh_(p.mesh()),
    p_(p),
    U_(U),
    T_(T),
    thermo_(thermo),
    // MRF_(MRF),
    rhoFlux_
    (
        IOobject
        (
            "phi",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (linearInterpolate(thermo_.rho()*U_) & mesh_.Sf())
    ),
    rhoUFlux_
    (
        IOobject
        (
            "rhoUFlux",
            mesh_.time().timeName(),
            mesh_,
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
            mesh_.time().timeName(),
            mesh_, 
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rhoFlux_*linearInterpolate(thermo.Cv()*T_ + 0.5*magSqr(U_))
    )//,
     
 
    // grad_( new newGradient(mesh_) ),

    // //     // Gradient of cell linear momentum
    // UGrad_(grad_->gradient(U_)),
    // pGrad_(grad_->gradient(p_)),
    // tGrad_(grad_->gradient(T))
    
    // meshPhi_
    // (
    //     IOobject
    //     (
    //         "meshPhi",
    //         mesh_.time().timeName(),
    //         mesh_,
    //         IOobject::NO_READ,
    //         IOobject::NO_WRITE
    //     ),
    //     mesh_,
    //     dimensionedScalar("0", dimVolume/dimTime, 0.0)
    // )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::numericFlux::computeFlux()
{
    // Get face-to-cell addressing: face area point from owner to neighbour
    const auto& owner = mesh_.owner();
    const auto& neighbour = mesh_.neighbour();

    // Get the face area vector
    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceScalarField& magSf = mesh_.magSf();
    
    // ALE mesh velocity + velocity due to MRF
        // Info << "evaluate meshPhi is evaluated"<<endl;
    surfaceScalarField mshPhi( meshPhi() ); 
    // Info << "meshPhi is evaluated"<<endl;
    // MRF_.makeAbsolute(mshPhi);

    // Thermodynamics
    const volScalarField Cv = thermo_.Cv();
    const volScalarField R  = thermo_.Cp() - Cv;

    surfaceScalarField pos_(IOobject("pos", mesh_), mesh_, dimensionedScalar("one", dimless, 1.0));
    surfaceScalarField neg_(IOobject("neg", mesh_), mesh_, dimensionedScalar("minusOne", dimless, -1.0));
   
    surfaceScalarField p_pos( fvc::interpolate(p_, pos_, "reconstruct(p)") );
    surfaceScalarField p_neg( fvc::interpolate(p_, neg_, "reconstruct(p)") );

    surfaceVectorField U_pos( fvc::interpolate(U_, pos_, "reconstruct(U)") );
    surfaceVectorField U_neg( fvc::interpolate(U_, neg_, "reconstruct(U)") );

    surfaceScalarField T_pos( fvc::interpolate(T_, pos_, "reconstruct(T)") );
    surfaceScalarField T_neg( fvc::interpolate(T_, neg_, "reconstruct(T)") );


    // surfaceScalarField p_pos
    // (
    //     IOobject
    //     (
    //         "p_pos",
    //         mesh_.time().timeName(),
    //         mesh_,
    //         IOobject::NO_READ,
    //         IOobject::NO_WRITE
    //     ),
    //     mesh_,
    //     dimensionedScalar("zero", dimPressure, 0.0)   // initialized to zero
    // );

    // surfaceScalarField p_neg
    // (
    //     IOobject
    //     (
    //         "p_neg",
    //         mesh_.time().timeName(),
    //         mesh_,
    //         IOobject::NO_READ,
    //         IOobject::NO_WRITE
    //     ),
    //     mesh_,
    //     dimensionedScalar("zero", dimPressure, 0.0)
    // );

    // surfaceVectorField U_pos
    // (
    //     IOobject
    //     (
    //         "U_pos",
    //         mesh_.time().timeName(),
    //         mesh_,
    //         IOobject::NO_READ,
    //         IOobject::NO_WRITE
    //     ),
    //     mesh_,
    //     dimensionedVector("zero", dimVelocity, vector::zero)
    // );

    // surfaceVectorField U_neg
    // (
    //     IOobject
    //     (
    //         "U_neg",
    //         mesh_.time().timeName(),
    //         mesh_,
    //         IOobject::NO_READ,
    //         IOobject::NO_WRITE
    //     ),
    //     mesh_,
    //     dimensionedVector("zero", dimVelocity, vector::zero)
    // );

    // surfaceScalarField T_pos
    // (
    //     IOobject
    //     (
    //         "T_pos",
    //         mesh_.time().timeName(),
    //         mesh_,
    //         IOobject::NO_READ,
    //         IOobject::NO_WRITE
    //     ),
    //     mesh_,
    //     dimensionedScalar("zero", dimTemperature, 0.0)
    // );

    // surfaceScalarField T_neg
    // (
    //     IOobject
    //     (
    //         "T_neg",
    //         mesh_.time().timeName(),
    //         mesh_,
    //         IOobject::NO_READ,
    //         IOobject::NO_WRITE
    //     ),
    //     mesh_,
    //     dimensionedScalar("zero", dimTemperature, 0.0)
    // );
    // Info<< "grad fields are created"<<endl;

    // UGrad_ = grad_->gradient(U_);
    // pGrad_ = grad_->gradient(p_);
    // tGrad_ = grad_->gradient(T_);
    // Info<< "grad_ is done"<<endl;
    //     // Reconstruction
    // grad_->reconstruct(U_, UGrad_, U_pos, U_neg);
    // grad_->reconstruct(p_, pGrad_, p_pos, p_neg);
    // grad_->reconstruct(T_, tGrad_, T_pos, T_neg);

    // Info<< "reconstruction is done"<<endl;
    // Calculate fluxes at internal faces
    forAll (owner, faceI)
    {
        const label own = owner[faceI];
        const label nei = neighbour[faceI];

        // Info<< "loop over fluxes "<<endl;
        // calculate fluxes with reconstructed primitive variables at faces
	    pFlux_ -> evaluateFlux
        (
            rhoFlux_[faceI],
            rhoUFlux_[faceI],
            rhoEFlux_[faceI],
            p_pos[faceI],  p_neg[faceI],
            U_pos[faceI],  U_neg[faceI],
            T_pos[faceI],  T_neg[faceI],
            R[own],        R[nei],
            Cv[own],       Cv[nei],
            Sf[faceI],
            magSf[faceI],
	        mshPhi[faceI]
        );
    }
    // Info<< "flux pFlux_ -> evaluateFlux  done"<<endl;

    // Update boundary field and values
    forAll (rhoFlux_.boundaryField(), patchi)
    {
        const fvPatch& curPatch = p_.boundaryField()[patchi].patch();

        // Fluxes
        fvsPatchScalarField& pRhoFlux  = rhoFlux_.boundaryField()[patchi];
        fvsPatchVectorField& pRhoUFlux = rhoUFlux_.boundaryField()[patchi];
        fvsPatchScalarField& pRhoEFlux = rhoEFlux_.boundaryField()[patchi];

        const scalarField& pCv = Cv.boundaryField()[patchi];
        const scalarField& pR  = R.boundaryField()[patchi];

        // Face areas
        const fvsPatchVectorField& pSf = Sf.boundaryField()[patchi];
        const fvsPatchScalarField& pMagSf = magSf.boundaryField()[patchi];
        const fvsPatchScalarField& pMshPhi = mshPhi.boundaryField()[patchi];

        if (curPatch.coupled())
        {
            // Patch fields
            const fvsPatchScalarField& pp_pos = p_pos.boundaryField()[patchi];
            const fvsPatchVectorField& pU_pos = U_pos.boundaryField()[patchi];
            const fvsPatchScalarField& pT_pos = T_pos.boundaryField()[patchi];
            
            const fvsPatchScalarField& pp_neg = p_neg.boundaryField()[patchi];
            const fvsPatchVectorField& pU_neg = U_neg.boundaryField()[patchi];
            const fvsPatchScalarField& pT_neg = T_neg.boundaryField()[patchi];
            
            forAll (curPatch, facei)
            {
                pFlux_ -> evaluateFlux
                (
                    pRhoFlux[facei],
                    pRhoUFlux[facei],
                    pRhoEFlux[facei],

                    pp_pos[facei],  pp_neg[facei],
                    pU_pos[facei],  pU_neg[facei],
                    pT_pos[facei],  pT_neg[facei],

                    pR[facei],  pR[facei],
                    pCv[facei], pCv[facei],
                    pSf[facei],
                    pMagSf[facei],
                    pMshPhi[facei]
                );
            }
        }
        else
        {
            const fvPatchScalarField& pp = p_.boundaryField()[patchi];
            const vectorField& pU = U_.boundaryField()[patchi];
            const scalarField& pT = T_.boundaryField()[patchi];

            forAll (pp, facei)
            {
                // Calculate fluxes
                pFlux_ -> evaluateFlux
                (
                    pRhoFlux[facei],
                    pRhoUFlux[facei],
                    pRhoEFlux[facei],
                    pp[facei],  pp[facei],
                    pU[facei],  pU[facei],
                    pT[facei],  pT[facei],
                    pR[facei],  pR[facei],
                    pCv[facei], pCv[facei],
                    pSf[facei],
                    pMagSf[facei],
		            pMshPhi[facei]
                );
            }
        }
    }

        // Info<< "flux boundary  done"<<endl;

}


Foam::tmp<Foam::surfaceScalarField> numericFlux::meshPhi() const

{
    if (this->mesh().moving()) 
    {
        // Info << "mesh is moving"<<endl;
        return  fvc::meshPhi(U_);
    } 

    return tmp<surfaceScalarField>
        (
            new surfaceScalarField
            (
                IOobject
                (
                "meshPhi",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
                ),
                this->mesh(),
                dimensionedScalar("0", dimVolume/dimTime, 0.0)
            )
        );
    
}


// ************************************************************************* //
