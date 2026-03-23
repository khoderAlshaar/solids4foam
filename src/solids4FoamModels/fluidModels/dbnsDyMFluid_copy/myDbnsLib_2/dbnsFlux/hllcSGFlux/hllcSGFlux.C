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

#include "hllcSGFlux.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::hllcSGFlux::evaluateFlux
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
    const scalar& magSf
) const
{

    // Step 1: decode left and right:
    // normal vector
    const vector normalVector = Sf/magSf;

    // Compute conservative variables assuming stiffened gas law

    // Density
    const scalar rhoLeft = (pLeft + pinf)/(Cv*(gamma-1)*TLeft);
    const scalar rhoRight = (pRight + pinf)/(Cv*(gamma-1)*TRight);

    // DensityVelocity
    const vector rhoULeft = rhoLeft*ULeft;
    const vector rhoURight = rhoRight*URight;

    // DensityTotalEnergy
    const scalar rhoELeft = rhoLeft*Cv*TLeft + pinf + rhoLeft*q + 0.5*rhoLeft*magSqr(ULeft);
    const scalar rhoERight = rhoRight*Cv*TRight + pinf + rhoRight*q + 0.5*rhoRight*magSqr(URight);

    // Compute left and right total enthalpies:
    const scalar HLeft = (rhoELeft + pLeft)/rhoLeft;
    const scalar HRight = (rhoERight + pRight)/rhoRight;

    // Compute UnLeft and UnRight (q_{l,r} = U_{l,r} \bullet n)
    const scalar UnLeft = (ULeft & normalVector);
    const scalar UnRight = (URight & normalVector);

    // Speed of sound, for left and right side, assuming stiffened gas
    const scalar aLeft =
        Foam::sqrt(max(0.0,gamma*(pLeft + pinf)/rhoLeft));

    const scalar aRight =
        Foam::sqrt(max(0.0,gamma*(pRight + pinf)/rhoRight));


    // Step 2:
    // needs rho_{l,r}, U_{l,r}, H_{l,r}, kappa_{l,r}, Gamma_{l,r}, q_{l,r}

    // Compute Roe weights
    const scalar rhoLeftSqrt = Foam::sqrt(max(0.0,rhoLeft));
    const scalar rhoRightSqrt = Foam::sqrt(max(0.0,rhoRight));

    const scalar wLeft = rhoLeftSqrt
        /stabilise((rhoLeftSqrt + rhoRightSqrt),VSMALL);

    const scalar wRight = 1 - wLeft;

    // Roe averaged velocity
    const vector UTilde = wLeft*ULeft + wRight*URight;

    // Roe averaged contravariant velocity
    const scalar contrUTilde = (UTilde & normalVector);

    // Roe averaged total enthalpy
    const scalar HTilde = wLeft*HLeft + wRight*HRight;

    // Speed of sound with Roe reconstruction values
    // TODO: not sure if the correct (flow speed) and kappa is used here
    const scalar aTilde =
        Foam::sqrt(max(0 ,(gamma - 1)*(HTilde - 0.5*magSqr(UTilde) - q)));

    // Step 3: compute signal speeds for face:
    const scalar SLeft  = min(UnLeft-aLeft, contrUTilde-aTilde);
    const scalar SRight = max(contrUTilde + aTilde, UnRight+aRight);

    const scalar SStar = (rhoRight*UnRight*(SRight-UnRight)
      - rhoLeft*UnLeft*(SLeft - UnLeft) + pLeft - pRight )/
        stabilise((rhoRight*(SRight-UnRight)-rhoLeft*(SLeft-UnLeft)),VSMALL);

    // Compute pressure in star region from the right side
    const scalar pStarRight =
        rhoRight*(UnRight - SRight)*(UnRight - SStar) + pRight;

    // Should be equal to the left side
    const scalar pStarLeft  =
        rhoLeft*(UnLeft -  SLeft)*(UnLeft - SStar) + pLeft;

    // Give a warning if this is not the case
    if (mag(pStarRight - pStarLeft) > 1e-6)
    {
        Info << "mag(pStarRight-pStarLeft) > VSMALL " << endl;
    }

    // Use pStarRight for pStar, as in theory, pStarRight == pStarLeft
    const scalar pStar = pStarRight;

    // Step 4: upwinding - compute states:
    scalar convectionSpeed = 0.0;
    scalar rhoState = 0.0;
    vector rhoUState = vector::zero;
    scalar rhoEState = 0.0;
    scalar pState = 0.0;

    if (pos(SLeft))
    {
        // compute F_l
        convectionSpeed = UnLeft;
        rhoState  = rhoLeft;
        rhoUState = rhoULeft;
        rhoEState = rhoELeft;
        pState = pLeft;
    }
    else if (pos(SStar))
    {
        scalar omegaLeft = scalar(1.0)/stabilise((SLeft - SStar), VSMALL);

        // Compute left star region
        convectionSpeed = SStar;
        rhoState  = omegaLeft*(SLeft - UnLeft)*rhoLeft;
        rhoUState = omegaLeft*((SLeft - UnLeft)*rhoULeft
        + (pStar - pLeft)*normalVector);
        rhoEState = omegaLeft*((SLeft - UnLeft)*rhoELeft
        - pLeft*UnLeft + pStar*SStar);
        pState = pStar;
    }
    else if (pos(SRight))
    {
        scalar omegaRight = scalar(1.0)/stabilise((SRight - SStar), VSMALL);

        // compute right star region
        convectionSpeed = SStar;
        rhoState  = omegaRight*(SRight - UnRight)*rhoRight;
        rhoUState = omegaRight*((SRight - UnRight)*rhoURight
        + (pStar - pRight)*normalVector);
        rhoEState = omegaRight*((SRight - UnRight)*rhoERight
        - pRight*UnRight + pStar*SStar);
        pState = pStar;
    }
    else if (neg(SRight))
    {
        // compute F_r
        convectionSpeed = UnRight;
        rhoState  = rhoRight;
        rhoUState = rhoURight;
        rhoEState = rhoERight;
        pState = pRight;
    }
    else
    {
        Info << "Error in HLLC Riemann solver" << endl;
    }

    rhoFlux  = (convectionSpeed*rhoState)*magSf;
    rhoUFlux = (convectionSpeed*rhoUState+pState*normalVector)*magSf;
    rhoEFlux = (convectionSpeed*(rhoEState+pState))*magSf;
}


// ************************************************************************* //
