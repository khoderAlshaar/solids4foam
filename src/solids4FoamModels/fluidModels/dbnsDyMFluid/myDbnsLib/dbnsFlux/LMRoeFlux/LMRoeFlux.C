/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

#include "LMRoeFlux.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(LMRoeFlux, 0);
    addToRunTimeSelectionTable(dbnsFlux, LMRoeFlux, dictionary);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::LMRoeFlux::evaluateFlux
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
    const scalar& RLeft,
    const scalar& RRight,
    const scalar& CvLeft,
    const scalar& CvRight,
    const vector& Sf,
    const scalar& magSf,
    const scalar& meshPhi,
    const tensor& R,
    const tensor& RTranspos,
    const scalar& fp1Left,
    const scalar& fp1Right
) const
{
//   if (mag(meshPhi)>0.0) 
//     {
//       FatalError
//         << "This dbnsFlux is not ready to run with moving meshes." << nl
//         << exit(FatalError);
//     };
    // cell face *normal* velocity w_n
    // const scalar w_n = 0.0;
   
    const scalar w_n = meshPhi / (magSf + VSMALL);
    // Step 1: decode rho left and right:
    scalar rhoLeft = pLeft/(RLeft*TLeft);
    scalar rhoRight = pRight/(RRight*TRight);

    // Decode left and right total energy:
    scalar eLeft = CvLeft*TLeft + 0.5*magSqr(ULeft);
    scalar eRight = CvRight*TRight + 0.5*magSqr(URight);

    // Adiabatic exponent is constant for ideal gas but if Cp=Cp(T)
    // it must be computed for each cell and evaluated at each face
    // through reconstruction
    const scalar kappaLeft = (CvLeft + RLeft)/CvLeft;
    const scalar kappaRight = (CvRight + RRight)/CvRight;

    // normal vector
    vector normalVector = Sf/magSf;

    // Compute left and right contravariant velocities:
    const scalar contrVLeft  = (ULeft & normalVector);
    const scalar contrVRight = (URight & normalVector);

    // Compute left and right total enthalpies:
    const scalar hLeft = eLeft + pLeft/rhoLeft;
    const scalar hRight = eRight + pRight/rhoRight;

    // Step 2: compute Roe averged quantities for face:
    const scalar rhoTilde = sqrt(max(rhoLeft*rhoRight, SMALL));

    // Some temporary variables:
    const scalar rhoLeftSqrt = sqrt(max(rhoLeft, SMALL));
    const scalar rhoRightSqrt = sqrt(max(rhoRight, SMALL));

    const scalar wLeft = rhoLeftSqrt/(rhoLeftSqrt + rhoRightSqrt);
    const scalar wRight = 1 - wLeft;

    const vector UTilde = ULeft*wLeft + URight*wRight;
    const scalar hTilde = hLeft*wLeft + hRight*wRight;
    const scalar qTildeSquare = magSqr(UTilde);
    const scalar kappaTilde = kappaLeft*wLeft + kappaRight*wRight;

    // Speed of sound
    const scalar cTilde =
        sqrt(max((kappaTilde - 1)*(hTilde - 0.5*qTildeSquare), SMALL));

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
    
    // speeds of sound (left and right states)
    const scalar cLeft = sqrt
    (
        max
        (
            (kappaLeft - 1)*(hLeft - 0.5*magSqr(ULeft)),
            SMALL
        )
    );

    const scalar cRight = sqrt
    (
        max
        (
            (kappaRight - 1)*(hRight - 0.5*magSqr(URight)),
            SMALL
        )
    );
    // // Local Mach numbers at left and right states (based on total velocity magnitude)
    const scalar ML = mag(ULeft) / (cLeft + VSMALL);
    const scalar MR = mag(URight) / (cRight + VSMALL);
    
    // Scaling factor: min(1, max(ML, MR)) as per Equation (6) in the paper
    const scalar Ma_local = max(ML, MR); 

    const scalar zeta = min(1.0, Ma_local) ;

    // -------------------------
    // Scaled velocity jumps for L2-Roe
    // -------------------------
    
    // Scale normal velocity jump (LMRoe + L2Roe)
    const scalar deltaContrV_star = zeta * deltaContrV;
    
    // Tangential velocity components
    const vector deltaU_nVec = deltaContrV * normalVector;
    const vector deltaU_t = deltaU - deltaU_nVec;
    
    // Scale tangential velocity jump (L2Roe enhancement)
    // const vector deltaU_t_star = zeta * deltaU_t;
    const vector deltaU_t_star = deltaU_t;

    // Roe and Pike - formulation
    const scalar r1 = (deltaP - rhoTilde * cTilde * deltaContrV_star) / (2.0*cTilde*cTilde);
    const scalar r3 = (deltaP + rhoTilde * cTilde * deltaContrV_star) / (2.0*cTilde*cTilde);
    
    // Entropy wave strength (unscaled)
    const scalar r2 = deltaRho - deltaP/(cTilde*cTilde);


    // Step 5: compute l vectors

    // rho row:
    const scalar l1rho = 1;
    const scalar l2rho = 1;
    const scalar l3rho = 0;
    const scalar l4rho = 1;

    // first U column
    const vector l1U = UTilde - cTilde*normalVector;

    // second U column
    const vector l2U = UTilde;

    // third U column
    const vector l3U = deltaU_t_star;

    // fourth U column
    const vector l4U = UTilde + cTilde*normalVector;

    // E row
    const scalar l1e = hTilde - cTilde*contrVTilde;
    const scalar l2e = 0.5*qTildeSquare;
    const scalar l3e = (UTilde & deltaU_t_star);  // Consistent with scaled jump
    const scalar l4e = hTilde + cTilde*contrVTilde;

    // Step 6: compute eigenvalues

    // derived from algebra by hand, only for Euler equation usefull
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
    // const scalar cLeft = sqrt
    // (
    //     max
    //     (
    //         (kappaLeft - 1)*(hLeft - 0.5*magSqr(ULeft)),
    //         SMALL
    //     )
    // );

    // const scalar cRight = sqrt
    // (
    //     max
    //     (
    //         (kappaRight - 1)*(hRight - 0.5*magSqr(URight)),
    //         SMALL
    //     )
    // );

    // First eigenvalue: U - c
    // scalar eps = 2*max(0,(UR - cRight) - (UL - cLeft));
    // if (lambda1 < eps)
    // {
    //     lambda1 = (sqr(lambda1) + sqr(eps))/(2.0*eps);
    // }

    // // Second eigenvalue: U
    // eps = 2*max(0, UR - UL);
    // if (lambda2 < eps)
    // {
    //     lambda2 = (sqr(lambda2) + sqr(eps))/(2.0*eps);
    // }

    // // Third eigenvalue: U + c
    // eps = 2*max(0,(UR + cRight) - (UL + cLeft));
    // if (lambda3 < eps)
    // {
    //     lambda3 = (sqr(lambda3) + sqr(eps))/(2.0*eps);
    // }
//----------------L2 Roe correction-----------------------
    const scalar lambda1L = contrVLeft - cLeft  - w_n;
    const scalar lambda3L = contrVLeft+ cLeft  - w_n;

    // Right state eigenvalues
    const scalar lambda1R = contrVRight - cRight - w_n;
    const scalar lambda3R = contrVRight + cRight - w_n;

    const scalar deltaLambda1 = max(lambda1R-lambda1L,0.0);
    const scalar deltaLabmda3 = max(lambda3R-lambda3L,0.0);

    if (lambda1 < 2*deltaLambda1)
    {
        lambda1 = (sqr(lambda1)/(4*deltaLambda1)) + deltaLambda1;

    }
    if (lambda3 < 2*deltaLabmda3)
    {
        lambda3 = (sqr(lambda3)/(4*deltaLabmda3)) + deltaLabmda3;

    }

    // Step 8: Compute flux differences

    // Components of deltaF1
    const scalar diffF11 = lambda1*r1*l1rho;
    const vector diffF124 = lambda1*r1*l1U;
    const scalar diffF15 = lambda1*r1*l1e;

    // Components of deltaF2
    const scalar diffF21 = lambda2*(r2*l2rho + rhoTilde*l3rho);
    const vector diffF224 = lambda2*(r2*l2U + rhoTilde*l3U);
    const scalar diffF25 = lambda2*(r2*l2e + rhoTilde*l3e);

    // Components of deltaF3
    const scalar diffF31 = lambda3*r3*l4rho;
    const vector diffF324 = lambda3*r3*l4U;
    const scalar diffF35 = lambda3*r3*l4e;

    // Step 9: compute left and right fluxes

    // Left flux 5-vector
    const scalar fluxLeft11 = rhoLeft*contrVLeft;
    const vector fluxLeft124 = ULeft*fluxLeft11 + normalVector*pLeft;
    const scalar fluxLeft15 = hLeft*fluxLeft11;

    // Right flux 5-vector
    const scalar fluxRight11 = rhoRight*contrVRight;
    const vector fluxRight124 = URight*fluxRight11 + normalVector*pRight;
    const scalar fluxRight15 = hRight*fluxRight11;

    // Step 10: compute face flux 5-vector
    const scalar flux1 =
        0.5*(fluxLeft11 + fluxRight11 - (rhoLeft+rhoRight)*w_n  - (diffF11 + diffF21 + diffF31));

    const vector flux24 =
        0.5*(fluxLeft124 + fluxRight124 - (rhoLeft*ULeft+rhoRight*URight)*w_n - (diffF124 + diffF224 + diffF324));

    const scalar flux5 =
        0.5*(fluxLeft15 + fluxRight15 - (rhoLeft*eLeft+rhoRight*eRight)*w_n - (diffF15 + diffF25 + diffF35));

    // Compute private data
    rhoFlux  = flux1*magSf;
    rhoUFlux = flux24*magSf;
    rhoEFlux = flux5*magSf;
}

// ************************************************************************* //
