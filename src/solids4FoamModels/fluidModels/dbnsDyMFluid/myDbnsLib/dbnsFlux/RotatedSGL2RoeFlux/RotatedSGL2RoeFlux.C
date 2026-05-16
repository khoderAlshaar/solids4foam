/*---------------------------------------------------------------------------*\
  L2-Roe low-dissipation Roe flux for low Mach numbers
  Implementation based on Oßwald et al. (2016) and Rieper (2011).

  Key corrections:
  - Local Mach number based on max(ML, MR), not Roe-averaged velocity
  - Scaling applied to velocity jumps in wave strength computation
  - Shock switch properly applied to both normal and tangential components
  - Physical flux evaluation unchanged
\*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*\
  LMRoe - Low Mach number fix for Roe's approximate Riemann solver
  Implementation based on Rieper (2011) "A low-Mach number fix for Roe's 
  approximate Riemann solver", Journal of Computational Physics 230 (2011).

  Key differences from L2Roe:
  - Only the normal velocity jump is scaled by the local Mach number
  - Tangential velocity jumps remain UNSCALED
  - This is sufficient to fix the low Mach accuracy problem
  - L2Roe additionally scales tangential jumps to reduce dissipation
    at high wavenumbers (relevant for LES/DES)
\*---------------------------------------------------------------------------*/

#include "RotatedSGL2RoeFlux.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(RotatedSGL2RoeFlux, 0);
    addToRunTimeSelectionTable(dbnsFlux, RotatedSGL2RoeFlux, dictionary);
}

void Foam::RotatedSGL2RoeFlux::evaluateFlux
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
    const tensor& RT,
    const scalar& fp1Left,
    const scalar& fp1Right
) const
{
    // cell face *normal* velocity w_n
    const scalar w_n = meshPhi / (magSf + VSMALL);

    // //! Step 1: decode rho left and right:

    // Density
    const scalar rhoLeft = (pLeft + pinf)/(Cv*(gamma-1)*TLeft);
    const scalar rhoRight = (pRight + pinf)/(Cv*(gamma-1)*TRight);

    const scalar cLeft = Foam::sqrt(max((gamma*(pLeft + pinf))/rhoLeft,SMALL));
    const scalar cRight = Foam::sqrt(max((gamma*(pRight + pinf))/rhoRight,SMALL));

    vector UL_hat = R & ULeft;
    vector UR_hat = R & URight;

    scalar unL = UL_hat.x();
    scalar unR = UR_hat.x();

    // DensityTotalEnergy
    const scalar rhoELeft = rhoLeft*Cv*TLeft + pinf + rhoLeft*q + 0.5*rhoLeft*magSqr(unL);
    const scalar rhoERight = rhoRight*Cv*TRight + pinf + rhoRight*q + 0.5*rhoRight*magSqr(unR);

    // Compute left and right total enthalpies:
    const scalar HLeft = (rhoELeft + pLeft)/rhoLeft;
    const scalar HRight = (rhoERight + pRight)/rhoRight;

    // normal vector
    // vector normalVector = Sf/magSf;
    
    // Step 2: compute Roe averged quantities for face:

    // Some temporary variables:
    const scalar rhoLeftSqrt = sqrt(max(rhoLeft, SMALL));
    const scalar rhoRightSqrt = sqrt(max(rhoRight, SMALL));

    const scalar wLeft = rhoLeftSqrt/(rhoLeftSqrt + rhoRightSqrt);
    const scalar wRight = 1 - wLeft;

    const scalar rhoTilde = sqrt(max(rhoLeft*rhoRight, SMALL));

    const vector UTilde = UL_hat*wLeft + UR_hat*wRight;

    const scalar UTildeX =UTilde.x();
    const scalar UTildeY =UTilde.y();
    const scalar UTildeZ =UTilde.z();

    const scalar HTilde  = HLeft*wLeft + HRight*wRight;
    const scalar qTildeSquare = magSqr(UTildeX);

    // Roe Speed of sound
    const scalar cTilde =
        Foam::sqrt(max(0 ,(gamma - 1)*(HTilde - 0.5*qTildeSquare - q)));

    // Step 3: compute primitive differences:
    const scalar deltaP = pRight - pLeft;
    const scalar deltaRho = rhoRight - rhoLeft;
    const vector deltaU = UR_hat - UL_hat;
    // const scalar deltaContrV = (deltaU & normalVector);

    const scalar ML = mag(unL) / (cLeft + VSMALL);
    const scalar MR = mag(unR) / (cRight + VSMALL);
    
    // Scaling factor: min(1, max(ML, MR)) as per Equation (6) in the paper
    const scalar Ma_local = max(ML, MR); //!  increase

    // // Info << "Ma_local= "<<Ma_local<<endl;
    const scalar zeta = min(1.0, Ma_local) ;
    // -------------------------
    // Eigenvalues (wave speeds)
    // -------------------------
    
    scalar lambda1 = mag(UTildeX - cTilde - w_n);
    scalar lambda234 = mag(UTildeX - w_n);
    scalar lambda5 = mag(UTildeX + cTilde - w_n);

    
    // Step 7: vanLeer entropy correction
    const scalar lambda1L = unL - cLeft  - w_n;
    const scalar lambda5L = unL+ cLeft  - w_n;

    // Right state eigenvalues
    const scalar lambda1R = unR - cRight - w_n;
    const scalar lambda5R = unR + cRight - w_n;

    const scalar deltaLambda1 = max(lambda1R-lambda1L,0.0);
    const scalar deltaLabmda5 = max(lambda5R-lambda5L,0.0);

    if (lambda1 < 2*deltaLambda1)
    {
        lambda1 = (sqr(lambda1)/(4*deltaLambda1)) + deltaLambda1;

    }
    if (lambda5 < 2*deltaLabmda5)
    {
        lambda5 = (sqr(lambda5)/(4*deltaLabmda5)) + deltaLabmda5;

    }


    // Acoustic wave strengths with SCALED normal velocity jump
    const scalar r1 = (deltaP - rhoTilde * cTilde * deltaU.x()*zeta ) / (2.0*cTilde*cTilde);

    const scalar r2 = deltaRho - deltaP/(cTilde*cTilde);   
    
    const scalar r3 =  rhoTilde*deltaU.y()*zeta;
    
    const scalar r4 =  rhoTilde*deltaU.z()*zeta;  
    
    const scalar r5 = (deltaP + rhoTilde * cTilde * deltaU.x()*zeta ) / (2.0*cTilde*cTilde);
    
    const vector cTildeI = cTilde * vector(1,0,0);
 
    const scalar l1rho = 1.0;
    const scalar l2rho = 1.0;
    const scalar l3rho = 0.0;
    const scalar l4rho = 0.0;
    const scalar l5rho = 1.0;

    const vector l1U = UTilde - cTildeI;
    const vector l2U = UTilde;
    const vector l3U = vector(0,1,0);  // Shear wave uses SCALED tangential jump
    const vector l4U = vector(0,0,1);  // Shear wave uses SCALED tangential jump
    const vector l5U = UTilde + cTildeI;

    const scalar l1e = HTilde  - cTilde*UTildeX;
    const scalar l2e = 0.5*qTildeSquare;
    // const scalar l2e = 0.5*(UTildeX*UTildeX);
    const scalar l3e = UTildeY;  // Consistent with scaled jump
    const scalar l4e = UTildeZ;  // Consistent with scaled jump
    const scalar l5e = HTilde  + cTilde*UTildeX;


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

    const scalar diffF21 = lambda234*r2*l2rho;
    const vector diffF224 = lambda234*r2*l2U;
    const scalar diffF25 = lambda234*r2*l2e;

    const scalar diffF31 = lambda234*r3*l3rho;
    const vector diffF324 = lambda234*r3*l3U;
    const scalar diffF35 = lambda234*r3*l3e;


    const scalar diffF41 = lambda234*r4*l4rho;
    const vector diffF424 = lambda234*r4*l4U;
    const scalar diffF45 = lambda234*r4*l4e;


    const scalar diffF51 = lambda5*r5*l5rho;
    const vector diffF524 = lambda5*r5*l5U;
    const scalar diffF55 = lambda5*r5*l5e;


  
    // Step 9: compute left and right fluxes
    // -------------------------
    // Physical fluxes (unchanged from standard Roe)
    // -------------------------
    
    const scalar fluxLeft11 = rhoLeft*unL;
    const vector fluxLeft124 = UL_hat*fluxLeft11 + pLeft*vector(1,0,0); //!UL_hat?
    const scalar fluxLeft15 = HLeft*fluxLeft11;
    // const scalar fluxLeft15 = (rhoELeft + pLeft)*contrVLeft;

    const scalar fluxRight11 = rhoRight*unR;
    const vector fluxRight124 = UR_hat*fluxRight11 + pRight*vector(1,0,0);//!UR_hat?
    const scalar fluxRight15 = HRight*fluxRight11;
    // const scalar fluxRight15 = (rhoERight + pRight)*contrVRight;

    // Step 10: compute face flux 5-vector
    // -------------------------
    // Face flux assembly (Roe flux with ALE correction)
    // -------------------------
    
    const scalar flux1 = 0.5*(fluxLeft11 + fluxRight11 
                              - (rhoLeft + rhoRight)*w_n 
                              - (diffF11 + diffF21 + diffF31+ diffF41+ diffF51));

    const vector flux24 = 0.5*(fluxLeft124 + fluxRight124 
                               - (rhoLeft*ULeft + rhoRight*URight)*w_n 
                               - (diffF124 + diffF224 + diffF324+ diffF424+ diffF524));

    const scalar flux5 = 0.5*(fluxLeft15 + fluxRight15 
                              - (rhoELeft + rhoERight )*w_n 
                              - (diffF15 + diffF25 + diffF35+ diffF45+ diffF55));

    // Scale by face area
    rhoFlux  = flux1*magSf;
    rhoUFlux = (RT & flux24)*magSf;
    rhoEFlux = flux5*magSf;
}

// ************************************************************************* //

// void Foam::RotatedSGL2RoeFlux::evaluateFlux_Wall
// (
//     scalar& rhoFlux,
//     vector& rhoUFlux,
//     scalar& rhoEFlux,
//     const scalar& pLeft,
//     const scalar& pRight,
//     const vector& ULeft,
//     const vector& URight,
//     const scalar& TLeft,
//     const scalar& TRight,
//     const scalar& RLeft,
//     const scalar& RRight,
//     const scalar& CvLeft,
//     const scalar& CvRight,
//     const vector& Sf,
//     const scalar& magSf,
//     const scalar& meshPhi,
//     const tensor& R,
//     const tensor& RT,
//     const scalar& fp1Left,
//     const scalar& fp1Right
// ) const
// {
//     // cell face *normal* velocity w_n
//     const scalar w_n = meshPhi / (magSf + VSMALL);

//     // //! Step 1: decode rho left and right:

//     // Density
//     const scalar rhoLeft = (pLeft + pinf)/(Cv*(gamma-1)*TLeft);
//     const scalar rhoRight = (pRight + pinf)/(Cv*(gamma-1)*TRight);

//     const scalar cLeft = Foam::sqrt(max((gamma*(pLeft + pinf))/rhoLeft,SMALL));
//     const scalar cRight = Foam::sqrt(max((gamma*(pRight + pinf))/rhoRight,SMALL));

//     vector UL_hat = R & ULeft;
//     vector UR_hat = R & URight;

//      UR_hat.x() = - UL_hat.x() ;

//     // DensityTotalEnergy
//     const scalar rhoELeft = rhoLeft*Cv*TLeft + pinf + rhoLeft*q + 0.5*rhoLeft*magSqr(UL_hat);
//     const scalar rhoERight = rhoRight*Cv*TRight + pinf + rhoRight*q + 0.5*rhoRight*magSqr(UR_hat);

//     // Compute left and right total enthalpies:
//     const scalar HLeft = (rhoELeft + pLeft)/rhoLeft;
//     const scalar HRight = (rhoERight + pRight)/rhoRight;

//     // normal vector
//     vector normalVector = Sf/magSf;
    
//     // Step 2: compute Roe averged quantities for face:

//     // Some temporary variables:
//     const scalar rhoLeftSqrt = sqrt(max(rhoLeft, SMALL));
//     const scalar rhoRightSqrt = sqrt(max(rhoRight, SMALL));

//     const scalar wLeft = rhoLeftSqrt/(rhoLeftSqrt + rhoRightSqrt);
//     const scalar wRight = 1 - wLeft;

//     const scalar rhoTilde = sqrt(max(rhoLeft*rhoRight, SMALL));

//     const vector UTilde = UL_hat*wLeft + UR_hat*wRight;

//     const scalar UTildeX =UTilde.x();
//     const scalar UTildeY =UTilde.y();
//     const scalar UTildeZ =UTilde.z();

//     const scalar HTilde  = HLeft*wLeft + HRight*wRight;
//     const scalar qTildeSquare = magSqr(UTilde);

//     // Roe Speed of sound
//     const scalar cTilde =
//         sqrt(max((gamma - 1)*(HTilde  - 0.5*qTildeSquare) - q , SMALL));

// //!! low mach number correction

//     const scalar MachLeft  = mag(ULeft) / (cLeft  + VSMALL);
//     const scalar MachRight = mag( URight) / (cRight + VSMALL);

// 	// scalar zeta = 1;
// 	scalar zeta = min(max(MachLeft,MachRight),1.0);

//     // scalar fp = pow3(min( (pLeft/pRight) , (pRight/pLeft)));
//     // scalar fp = min(fp1Left,fp1Right);
    
//     //  scalar zeta = 1- (1-zeta1) *fp ;

// //    if (mesh_.time().outputTime())
// //     {
// //     Info << "zeta1= "<<zeta1 <<",   zeta= " << zeta <<", fp= "<< fp <<endl;
// //     }
//     // zeta = 0.1;//1- (1-zeta) *fp ;

//     // Step 1 — Left and Right Acoustic Eigenvalues
//     // Left state eigenvalues
//     const scalar lambda1L = UL_hat.x() - cLeft  - w_n;
//     const scalar lambda5L = UL_hat.x() + cLeft  - w_n;

//     // Right state eigenvalues
//     const scalar lambda1R = UR_hat.x() - cRight - w_n;
//     const scalar lambda5R = UR_hat.x() + cRight - w_n;

//     const scalar deltaLambda1 = max(lambda1R-lambda1L,0.0);
//     const scalar deltaLabmda5 = max(lambda5R-lambda5L,0.0);


//     // Step 3: compute primitive differences:
//     const scalar deltaP = pRight - pLeft;
//     const scalar deltaRho = rhoRight - rhoLeft;
//     const vector deltaU = UR_hat - UL_hat;
//     // const scalar deltaContrV = (deltaU & normalVector);


//     // -------------------------
//     // Eigenvalues (wave speeds)
//     // -------------------------
    
//     scalar lambda1 = mag(UTildeX - cTilde - w_n);
//     scalar lambda234 = mag(UTildeX - w_n);
//     scalar lambda5 = mag(UTildeX + cTilde - w_n);

    
//     // Step 7: vanLeer entropy correction
//     if (lambda1 < 2*deltaLambda1)
//     {
//         lambda1 = (sqr(lambda1)/(4*deltaLambda1)) + deltaLambda1;

//     }
//     if (lambda5 < 2*deltaLabmda5)
//     {
//         lambda5 = (sqr(lambda5)/(4*deltaLabmda5)) + deltaLabmda5;

//     }


//     // Acoustic wave strengths with SCALED normal velocity jump
//     const scalar r1 = (deltaP - rhoTilde * cTilde * deltaU.x()*zeta ) / (2.0*cTilde*cTilde);

//     const scalar r2 = deltaRho - deltaP/(cTilde*cTilde);   
    
//     const scalar r3 =  rhoTilde*deltaU.y()*zeta;
    
//     const scalar r4 =  rhoTilde*deltaU.z()*zeta;  
    
//     const scalar r5 = (deltaP + rhoTilde * cTilde * deltaU.x()*zeta ) / (2.0*cTilde*cTilde);
    
//     const vector cTildeI = cTilde * vector(1,0,0);
 
//     const scalar l1rho = 1.0;
//     const scalar l2rho = 1.0;
//     const scalar l3rho = 0.0;
//     const scalar l4rho = 0.0;
//     const scalar l5rho = 1.0;

//     const vector l1U = UTilde - cTildeI;
//     const vector l2U = UTilde;
//     const vector l3U = vector(R.yx(), R.yy(), R.yz()); //vector(0,1,0);  // Shear wave uses SCALED tangential jump
//     const vector l4U = vector(R.zx(), R.zy(), R.zz()); //vector(0,0,1);  // Shear wave uses SCALED tangential jump
//     const vector l5U = UTilde + cTildeI;

//     const scalar l1e = HTilde  - cTilde*UTildeX;
//     const scalar l2e = 0.5*qTildeSquare;
//     // const scalar l2e = 0.5*(UTildeX*UTildeX);
//     const scalar l3e = UTildeY;  // Consistent with scaled jump
//     const scalar l4e = UTildeZ;  // Consistent with scaled jump
//     const scalar l5e = HTilde  + cTilde*UTildeX;


//     // -------------------------
//     // Flux difference components
//     // -------------------------


//     // Step 8: Compute flux differences

//     // -------------------------
//     // Flux difference components
//     // -------------------------
    
//     const scalar diffF11 = lambda1*r1*l1rho;
//     const vector diffF124 = lambda1*r1*l1U;
//     const scalar diffF15 = lambda1*r1*l1e;

//     const scalar diffF21 = lambda234*r2*l2rho;
//     const vector diffF224 = lambda234*r2*l2U;
//     const scalar diffF25 = lambda234*r2*l2e;

//     const scalar diffF31 = lambda234*r3*l3rho;
//     const vector diffF324 = lambda234*r3*l3U;
//     const scalar diffF35 = lambda234*r3*l3e;


//     const scalar diffF41 = lambda234*r4*l4rho;
//     const vector diffF424 = lambda234*r4*l4U;
//     const scalar diffF45 = lambda234*r4*l4e;


//     const scalar diffF51 = lambda5*r5*l5rho;
//     const vector diffF524 = lambda5*r5*l5U;
//     const scalar diffF55 = lambda5*r5*l5e;


  
//     // Step 9: compute left and right fluxes
//     // -------------------------
//     // Physical fluxes (unchanged from standard Roe)
//     // -------------------------
    
//     const scalar fluxLeft11 = rhoLeft*UL_hat.x();
//     const vector fluxLeft124 = UL_hat*fluxLeft11 + pLeft*vector(1,0,0); //!normalVector*pRight;
//     const scalar fluxLeft15 = HLeft*fluxLeft11;
//     // const scalar fluxLeft15 = (rhoELeft + pLeft)*contrVLeft;

//     const scalar fluxRight11 = rhoRight*UR_hat.x();
//     const vector fluxRight124 = UR_hat*fluxRight11 + pRight*vector(1,0,0);
//     const scalar fluxRight15 = HRight*fluxRight11;
//     // const scalar fluxRight15 = (rhoERight + pRight)*contrVRight;

//     // Step 10: compute face flux 5-vector
//     // -------------------------
//     // Face flux assembly (Roe flux with ALE correction)
//     // -------------------------
    
//     const scalar flux1 = 0.5*(fluxLeft11 + fluxRight11 
//                               - (rhoLeft + rhoRight)*w_n 
//                               - (diffF11 + diffF21 + diffF31+ diffF41+ diffF51));

//     const vector flux24 = 0.5*(fluxLeft124 + fluxRight124 
//                                - (rhoLeft*ULeft + rhoRight*URight)*w_n 
//                                - (diffF124 + diffF224 + diffF324+ diffF424+ diffF524));

//     const scalar flux5 = 0.5*(fluxLeft15 + fluxRight15 
//                               - (rhoELeft + rhoERight )*w_n 
//                               - (diffF15 + diffF25 + diffF35+ diffF45+ diffF55));

//     // Scale by face area
//     rhoFlux  = flux1*magSf;
//     rhoUFlux = (RT & flux24)*magSf;
//     rhoEFlux = flux5*magSf;
// }

// ************************************************************************* //
