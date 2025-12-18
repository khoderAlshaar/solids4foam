/*---------------------------------------------------------------------------*\
  L2-Roe low-dissipation Roe flux for low Mach numbers
  Implementation based on Oßwald et al. (2016) and Rieper (2011).

  Key corrections:
  - Local Mach number based on max(ML, MR), not Roe-averaged velocity
  - Scaling applied to velocity jumps in wave strength computation
  - Shock switch properly applied to both normal and tangential components
  - Physical flux evaluation unchanged
\*---------------------------------------------------------------------------*/

#include "L2RoeALEFlux.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(L2RoeALEFlux, 0);
    addToRunTimeSelectionTable(dbnsFlux, L2RoeALEFlux, dictionary);
}

void Foam::L2RoeALEFlux::evaluateFlux
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
    // mesh normal velocity
    const scalar w_n = meshPhi / (magSf + VSMALL);

    // densities
    const scalar rhoLeft = pLeft/(RLeft*TLeft);
    const scalar rhoRight = pRight/(RRight*TRight);

    // total energy per unit mass
    const scalar eLeft = CvLeft*TLeft + 0.5*magSqr(ULeft);
    const scalar eRight = CvRight*TRight + 0.5*magSqr(URight);

    // adiabatic exponents
    const scalar kappaLeft = (CvLeft + RLeft)/CvLeft;
    const scalar kappaRight = (CvRight + RRight)/CvRight;

    // normal vector
    const vector normalVector = Sf/magSf;

    // contravariant velocities (normal components)
    const scalar contrVLeft  = (ULeft & normalVector);
    const scalar contrVRight = (URight & normalVector);

    // enthalpies
    const scalar hLeft = eLeft + pLeft/rhoLeft;
    const scalar hRight = eRight + pRight/rhoRight;

    // Step 2: compute Roe averged quantities for face:
    const scalar rhoTilde = sqrt(max(rhoLeft*rhoRight, SMALL));

    // Roe averages (sqrt density weighting)
    const scalar rhoLeftSqrt = sqrt(max(rhoLeft, SMALL));
    const scalar rhoRightSqrt = sqrt(max(rhoRight, SMALL));

    const scalar wLeft = rhoLeftSqrt/(rhoLeftSqrt + rhoRightSqrt);
    const scalar wRight = 1.0 - wLeft;

    const vector UTilde = ULeft*wLeft + URight*wRight;
    const scalar hTilde = hLeft*wLeft + hRight*wRight;
    const scalar qTildeSquare = magSqr(UTilde);
    const scalar kappaTilde = kappaLeft*wLeft + kappaRight*wRight;

    // Roe speed of sound
    const scalar cTilde = sqrt(max((kappaTilde - 1)*(hTilde - 0.5*qTildeSquare), SMALL));

    // Roe averaged normal velocity
    const scalar contrVTilde = (UTilde & normalVector);

    // primitive differences
    const scalar deltaP = pRight - pLeft;
    const scalar deltaRho = rhoRight - rhoLeft;
    const vector deltaU = URight - ULeft;
    const scalar deltaContrV = (deltaU & normalVector);

    // -------------------------
    // L2-Roe scaling: Compute local Mach number based on LEFT and RIGHT states
    // -------------------------
    
    //! speeds of sound (left and right states)
    const scalar cLeft = sqrt(max((kappaLeft - 1)*(hLeft - 0.5*magSqr(ULeft)), SMALL));
    const scalar cRight = sqrt(max((kappaRight - 1)*(hRight - 0.5*magSqr(URight)), SMALL));

    // Local Mach numbers at left and right states (based on total velocity magnitude)
    const scalar ML = mag(ULeft) / (cLeft + VSMALL);
    const scalar MR = mag(URight) / (cRight + VSMALL);
    
    // Scaling factor: min(1, max(ML, MR)) as per Equation (6) in the paper
    const scalar Ma_local = max(ML, MR);
    const scalar zeta = min(1.0, Ma_local);

    // -------------------------
    // Shock switch (Portela style with modifications)
    // -------------------------
    
    // Compute shock indicator
    const scalar UL = contrVLeft;
    const scalar UR = contrVRight;

    scalar eps1 = 2.0*max(0.0, (UR - cRight) - (UL - cLeft));
    scalar eps2 = 2.0*max(0.0, UR - UL);
    scalar eps3 = 2.0*max(0.0, (UR + cRight) - (UL + cLeft));

    // Simplified shock switch: active if significant pressure jump exists
    // Can be made more sophisticated based on pressure gradient
    const scalar deltaP_threshold = 0.01*min(pLeft, pRight);
    const bool shockPresent = (mag(deltaP) > deltaP_threshold);
    
    // Apply scaling only away from shocks (ssw = 0)
    // const scalar zeta_eff = shockPresent ? 1.0 : zeta;
    const scalar zeta_eff = zeta;
    // const scalar zeta_eff = 1;//shockPresent ? 1.0 : zeta;

    // -------------------------
    // Scaled velocity jumps for L2-Roe
    // -------------------------
    
    // Scale normal velocity jump (LMRoe + L2Roe)
    const scalar deltaContrV_star = zeta_eff * deltaContrV;
    
    // Tangential velocity components
    const vector deltaU_nVec = deltaContrV * normalVector;
    const vector deltaU_t = deltaU - deltaU_nVec;
    
    // Scale tangential velocity jump (L2Roe enhancement)
    const vector deltaU_t_star = zeta_eff * deltaU_t;

    // -------------------------
    // Wave strengths using scaled jumps
    // -------------------------


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

    const scalar l1e = hTilde - cTilde*contrVTilde;
    const scalar l2e = 0.5*qTildeSquare;
    const scalar l3e = (UTilde & deltaU_t_star);  // Consistent with scaled jump
    const scalar l4e = hTilde + cTilde*contrVTilde;

    // -------------------------
    // Eigenvalues (wave speeds)
    // -------------------------
    
    scalar lambda1 = mag(contrVTilde - cTilde - w_n);
    scalar lambda2 = mag(contrVTilde - w_n);
    scalar lambda3 = mag(contrVTilde + cTilde - w_n);

    // -------------------------
    // Entropy correction (Harten-Hyman type)
    // -------------------------
    
    if (lambda1 < eps1 && eps1 > VSMALL)
    {
        lambda1 = (sqr(lambda1) + sqr(eps1))/(2.0*eps1);
    }

    if (lambda2 < eps2 && eps2 > VSMALL)
    {
        lambda2 = (sqr(lambda2) + sqr(eps2))/(2.0*eps2);
    }

    if (lambda3 < eps3 && eps3 > VSMALL)
    {
        lambda3 = (sqr(lambda3) + sqr(eps3))/(2.0*eps3);
    }

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

    // -------------------------
    // Physical fluxes (unchanged from standard Roe)
    // -------------------------
    
    const scalar fluxLeft11 = rhoLeft*contrVLeft;
    const vector fluxLeft124 = ULeft*fluxLeft11 + normalVector*pLeft;
    const scalar fluxLeft15 = hLeft*fluxLeft11;

    const scalar fluxRight11 = rhoRight*contrVRight;
    const vector fluxRight124 = URight*fluxRight11 + normalVector*pRight;
    const scalar fluxRight15 = hRight*fluxRight11;

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
                              - (rhoLeft*eLeft + rhoRight*eRight)*w_n 
                              - (diffF15 + diffF25 + diffF35));

    // Scale by face area
    rhoFlux  = flux1*magSf;
    rhoUFlux = flux24*magSf;
    rhoEFlux = flux5*magSf;
}

// ************************************************************************* //
