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

#include "SGL2RoeALEFlux.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::SGL2RoeALEFlux::evaluateFlux
(
    scalar& rhoFlux,
    vector& rhoUFlux,
    scalar& rhoEFlux,
    const scalar& pLeft,
    const scalar& pRight,
    const vector& ULeft,
    const vector& URight,
    const scalar& TLeft,
    const scalar& TRight,
    const scalar& q,
    const scalar& pinf,
    const scalar& gamma,
    const scalar& Cv,
    const vector& Sf,
    const scalar& magSf,
    const scalar& meshPhi,
    const scalar& fp1Left,
    const scalar& fp1Right
) const
{

    // Step 1: decode left and right:
    // normal vector
    const vector normalVector = Sf/magSf;

    // cell face *normal* velocity w_n
    const scalar w_n = meshPhi / (magSf + VSMALL);

    // Compute conservative variables assuming stiffened gas law

    // Density
    const scalar rhoLeft = (pLeft + pinf)/(Cv*(gamma-1)*TLeft);
    const scalar rhoRight = (pRight + pinf)/(Cv*(gamma-1)*TRight);

    // DensityTotalEnergy
    const scalar rhoELeft = rhoLeft*Cv*TLeft + pinf + rhoLeft*q + 0.5*rhoLeft*magSqr(ULeft);
    const scalar rhoERight = rhoRight*Cv*TRight + pinf + rhoRight*q + 0.5*rhoRight*magSqr(URight);

    // Compute left and right total enthalpies:
    const scalar HLeft = (rhoELeft + pLeft)/rhoLeft;
    const scalar HRight = (rhoERight + pRight)/rhoRight;

    // Compute left and right contravariant velocities:
    const scalar contrVLeft  = (ULeft & normalVector);
    const scalar contrVRight = (URight & normalVector);

    // Step 2: compute Roe averged quantities for face:
    const scalar rhoTilde = sqrt(max(rhoLeft*rhoRight, SMALL));

    // Some temporary variables:
    const scalar rhoLeftSqrt = sqrt(max(rhoLeft, SMALL));
    const scalar rhoRightSqrt = sqrt(max(rhoRight, SMALL));

    const scalar wLeft = rhoLeftSqrt/(rhoLeftSqrt + rhoRightSqrt);
    const scalar wRight = 1 - wLeft;

    const vector UTilde = ULeft*wLeft + URight*wRight;
    const scalar HTilde  = HLeft*wLeft + HRight*wRight;
    const scalar qTildeSquare = magSqr(UTilde);
    // const scalar kappaTilde = kappaLeft*wLeft + kappaRight*wRight;

    //! Roe Speed of sound (-q)
    const scalar cTilde =
        Foam::sqrt(max(0 ,(gamma - 1)*(HTilde - 0.5*magSqr(UTilde) - q)));

    // Roe averaged contravariant velocity
    const scalar contrVTilde = (UTilde & normalVector);

    // Step 3: compute primitive differences:
    const scalar deltaP = pRight - pLeft;
    const scalar deltaRho = rhoRight - rhoLeft;
    const vector deltaU = URight - ULeft;
    const scalar deltaContrV = (deltaU & normalVector);

    // Step 4: compute wave strengths:

    // -------------------------
    // L2-Roe scaling: Compute local Mach number based on LEFT and RIGHT states
    // -------------------------
    
    // // speeds of sound (left and right states)
    const scalar cLeft = Foam::sqrt(max((gamma*(pLeft + pinf))/rhoLeft,SMALL));
    const scalar cRight = Foam::sqrt(max((gamma*(pRight + pinf))/rhoRight,SMALL));
    
    // // Local Mach numbers at left and right states (based on total velocity magnitude)
    const scalar ML = mag(ULeft) / (cLeft + VSMALL);
    const scalar MR = mag(URight) / (cRight + VSMALL);
    
    // Scaling factor: min(1, max(ML, MR)) as per Equation (6) in the paper
    const scalar Ma_local = max(ML, MR); //!  increase

    // // Info << "Ma_local= "<<Ma_local<<endl;
    const scalar zeta = min(1.0, Ma_local) ;

    // const vector UTilde_normal = contrVTilde * normalVector;
    // const vector UTilde_tangent = UTilde - UTilde_normal;
    // const scalar VTildeMag = mag(UTilde_tangent);

    // -------------------------
    // Scaled velocity jumps for L2-Roe
    // -------------------------
    
    // Scale normal velocity jump (LMRoe + L2Roe)
    const scalar deltaContrV_star = zeta * deltaContrV;
    
    // Tangential velocity components
    const vector deltaU_nVec = deltaContrV * normalVector;
    const vector deltaU_t = deltaU - deltaU_nVec;
    
    // Scale tangential velocity jump (L2Roe enhancement)
    const vector deltaU_t_star = zeta * deltaU_t;
    // const vector deltaU_t_star = deltaU_t;

        // -------------------------
    // Wave strengths using scaled jumps
    // -------------------------
    // const scalar rhoTilde = sqrt(max(rhoLeft*rhoRight, SMALL));

    // Acoustic wave strengths with SCALED normal velocity jump
    const scalar r1 = (deltaP - rhoTilde * cTilde * deltaContrV_star) / (2.0*cTilde*cTilde);
    const scalar r3 = (deltaP + rhoTilde * cTilde * deltaContrV_star) / (2.0*cTilde*cTilde);
    
    // Entropy wave strength (unscaled)
    const scalar r2 = deltaRho - deltaP/(cTilde*cTilde);


    // -------------------------
    // Eigenvectors (right eigenvectors of Roe matrix)
    // -------------------------
    
    const scalar l1rho = 1.0;
    const scalar l2rho = 1.0;
    const scalar l3rho = 0.0;
    const scalar l4rho = 1.0;

    const vector l1U = UTilde - cTilde*normalVector;
    const vector l2U = UTilde;
    const vector l3U = deltaU_t_star;  // Shear wave uses SCALED tangential jump
    const vector l4U = UTilde + cTilde*normalVector;

    const scalar l1e = HTilde  - cTilde*contrVTilde;
    const scalar l2e = 0.5*qTildeSquare;
    const scalar l3e = (UTilde & deltaU_t_star);  // Consistent with scaled jump
    const scalar l4e = HTilde  + cTilde*contrVTilde;

    // Step 6: compute eigenvalues

    // -------------------------
    // Eigenvalues (wave speeds)
    scalar lambda1 = mag(contrVTilde - cTilde - w_n);
    scalar lambda2 = mag(contrVTilde - w_n);
    scalar lambda3 = mag(contrVTilde + cTilde - w_n);

    // Step 7: check for Harten entropy correction

//     const scalar eps = 0.1*cTilde; //adjustable parameter

//     if (lambda1 < eps || lambda2 < eps || lambda3 < eps)
//     {
//         lambda1 = (sqr(lambda1) + sqr(eps))/(2.0*eps);
//         lambda2 = (sqr(lambda2) + sqr(eps))/(2.0*eps);
//         lambda3 = (sqr(lambda3) + sqr(eps))/(2.0*eps);
//     }

    // Step 7a: Alternative entropy correction: Felipe Portela, 9/Oct/2013

    const scalar UL = ULeft & normalVector;
    const scalar UR = URight & normalVector;

    // First eigenvalue: U - c
    scalar eps = 2*max(0,(UR - cRight - w_n) - (UL - cLeft - w_n));
    if (lambda1 < eps)
    {
        lambda1 = (sqr(lambda1) + sqr(eps))/(2.0*eps);
    }

    // Second eigenvalue: U
    eps = 2*max(0, UR - UL - w_n);
    if (lambda2 < eps)
    {
        lambda2 = (sqr(lambda2) + sqr(eps))/(2.0*eps);
    }

    // Third eigenvalue: U + c
    eps = 2*max(0,(UR + cRight - w_n) - (UL + cLeft- w_n));
    if (lambda3 < eps)
    {
        lambda3 = (sqr(lambda3) + sqr(eps))/(2.0*eps);
    }
    // -------------------------
    // Flux difference components
    // -------------------------


    // Step 8: Compute flux differences

    // -------------------------
    // Flux difference components
    // -------------------------
    
    const scalar diffF11 = lambda1*r1*l1rho;
    const vector diffF124 = lambda1*r1*l1U;
    const scalar diffF15 = lambda1*r1*l1e;

    const scalar diffF21 = lambda2*(r2*l2rho + rhoTilde*l3rho);
    const vector diffF224 = lambda2*(r2*l2U + rhoTilde*l3U);
    const scalar diffF25 = lambda2*(r2*l2e + rhoTilde*l3e);

    const scalar diffF31 = lambda3*r3*l4rho;
    const vector diffF324 = lambda3*r3*l4U;
    const scalar diffF35 = lambda3*r3*l4e;
  
    // Step 9: compute left and right fluxes
    // -------------------------
    // Physical fluxes (unchanged from standard Roe)
    // -------------------------
    
    const scalar fluxLeft11 = rhoLeft*contrVLeft;
    const vector fluxLeft124 = ULeft*fluxLeft11 + normalVector*pLeft;
    const scalar fluxLeft15 = HLeft*fluxLeft11;
    // const scalar fluxLeft15 = (rhoELeft + pLeft)*contrVLeft;

    const scalar fluxRight11 = rhoRight*contrVRight;
    const vector fluxRight124 = URight*fluxRight11 + normalVector*pRight;
    const scalar fluxRight15 = HRight*fluxRight11;
    // const scalar fluxRight15 = (rhoERight + pRight)*contrVRight;

    // Step 10: compute face flux 5-vector
    // -------------------------
    // Face flux assembly (Roe flux with ALE correction)
    // -------------------------
    
    const scalar flux1 = 0.5*(fluxLeft11 + fluxRight11 
                              - (rhoLeft + rhoRight)*w_n 
                              - (diffF11 + diffF21 + diffF31));

    const vector flux24 = 0.5*(fluxLeft124 + fluxRight124 
                               - (rhoLeft*ULeft + rhoRight*URight)*w_n 
                               - (diffF124 + diffF224 + diffF324));

    const scalar flux5 = 0.5*(fluxLeft15 + fluxRight15 
                              - (rhoELeft + rhoERight )*w_n 
                              - (diffF15 + diffF25 + diffF35));

    // Scale by face area
    rhoFlux  = flux1*magSf;
    rhoUFlux = flux24*magSf;
    rhoEFlux = flux5*magSf;
}

// ************************************************************************* //
