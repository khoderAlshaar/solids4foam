/*---------------------------------------------------------------------------*\
  L2-Roe low-dissipation Roe flux for low Mach numbers
  Implementation based on Oßwald et al. (2016) and Rieper (2011).

  Key corrections:
  - Local Mach number based on max(ML, MR), not Roe-averaged velocity
  - Scaling applied to velocity jumps in wave strength computation
  - Shock switch properly applied to both normal and tangential components
  - Physical flux evaluation unchanged
\*---------------------------------------------------------------------------*/

#include "roeALELowMaFlux.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(roeALELowMaFlux, 0);
    addToRunTimeSelectionTable(dbnsFlux, roeALELowMaFlux, dictionary);
}

void Foam::roeALELowMaFlux::evaluateFlux
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
    const scalar& meshPhi
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
    // -------------------------
    // L2-Roe scaling: Compute local Mach number based on LEFT and RIGHT states
    // -------------------------
    //     const scalar UL = contrVLeft;
    // const scalar UR = contrVRight;
    
    // speeds of sound (left and right states)
    const scalar cLeft = sqrt(max((kappaLeft - 1)*(hLeft - 0.5*magSqr(ULeft)), SMALL));
    const scalar cRight = sqrt(max((kappaRight - 1)*(hRight - 0.5*magSqr(URight)), SMALL));


    // Local Mach numbers at left and right states (based on total velocity magnitude)
    const scalar ML = mag(ULeft) / (cLeft + VSMALL);
    const scalar MR = mag(URight) / (cRight + VSMALL);
    
    // Scaling factor: min(1, max(ML, MR)) as per Equation (6) in the paper
    const scalar Ma_local = max(ML, MR);
    const scalar zeta = min(1.0, Ma_local);


    // const vector UTilde_normal = contrVTilde * normalVector;
    // const vector UTilde_tangent = UTilde - UTilde_normal;
    // const scalar VTildeMag = mag(UTilde_tangent);
    
    // // Local Mach number as per Rieper (2011), Eq. (3.16):
    // // Ma_local = (|U_n| + |V_t|) / a
    // const scalar Ma_local = (mag(UTilde_normal) + mag(UTilde_tangent)) / (cTilde + VSMALL);
    
    // // Scaling factor: min(Ma_local, 1)
    // const scalar zeta =  min(1.0, Ma_local);

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

    // -------------------------
    // Wave strengths using scaled jumps
    // -------------------------

    // Acoustic wave strengths with SCALED normal velocity jump
    const scalar r1 = (deltaP - rhoTilde * cTilde * deltaContrV_star) / (2.0*sqr(cTilde));
    const scalar r3 = (deltaP + rhoTilde * cTilde * deltaContrV_star) / (2.0*sqr(cTilde));
    
    // Entropy wave strength (unscaled)
    const scalar r2 = deltaRho - deltaP/(2.0*sqr(cTilde));

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

    const scalar l1e = hTilde - cTilde*contrVTilde;
    const scalar l2e = 0.5*qTildeSquare;
    const scalar l3e = (UTilde & deltaU_t_star);  // Consistent with scaled jump
    const scalar l4e = hTilde + cTilde*contrVTilde;

    // -------------------------
    // Eigenvalues (wave speeds)
    // -------------------------
 
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
    scalar eps = 2*max(0,(UR - cRight) - (UL - cLeft));
    if (lambda1 < eps)
    {
        lambda1 = (sqr(lambda1) + sqr(eps))/(2.0*eps);
    }

    // Second eigenvalue: U
    eps = 2*max(0, UR - UL);
    if (lambda2 < eps)
    {
        lambda2 = (sqr(lambda2) + sqr(eps))/(2.0*eps);
    }

    // Third eigenvalue: U + c
    eps = 2*max(0,(UR + cRight) - (UL + cLeft));
    if (lambda3 < eps)
    {
        lambda3 = (sqr(lambda3) + sqr(eps))/(2.0*eps);
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
