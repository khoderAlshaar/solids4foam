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

#include "solidRoeFlux.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(solidRoeFlux, 0);
    addToRunTimeSelectionTable(solidFlux, solidRoeFlux, dictionary);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::solidRoeFlux::evaluateFlux
(
    vector& lmFlux,
    tensor& FFlux,
    const tensor& S_lm,
    const tensor& S_t,
    const vector& t_M,
    const vector& t_P,
    const vector& lm_P,
    const vector& lm_M,
    const tensor& F_P,
    const tensor& F_M,
    const tensor& P_P,
    const tensor& P_M,
    const tensor& R,
    const tensor& RT,
    const scalar rho,
    const scalar lambda,
    const scalar Up,
    const scalar Us,
    const scalar magSf,
    const vector Sf,
    const vector N
) const
{
    // Rotate states
    const vector lm_Mh = R & lm_M;
    const vector lm_Ph = R & lm_P;

    const tensor F_Mh = (R & F_M) & RT;
    const tensor F_Ph = (R & F_P) & RT;

    const tensor P_Mh = (R & P_M) & RT;
    const tensor P_Ph = (R & P_P) & RT;

    // Differences
    const vector dLm = lm_Ph - lm_Mh;
    const tensor dF  = F_Ph  - F_Mh;

    // Central flux (hat system)
    vector fluxLmHat
    (
        0.5*(P_Ph.xx() + P_Mh.xx()),
        0.5*(P_Ph.yx() + P_Mh.yx()),
        0.5*(P_Ph.zx() + P_Mh.zx())
    );

    const scalar invRho = 1.0/rho;

    tensor fluxFHat(tensor::zero);
    fluxFHat.xx() = invRho * 0.5*(lm_Ph.x() + lm_Mh.x());
    fluxFHat.yx() = invRho * 0.5*(lm_Ph.y() + lm_Mh.y());
    fluxFHat.zx() = invRho * 0.5*(lm_Ph.z() + lm_Mh.z());

    // Cached sums
    const scalar trF  = dF.yy() + dF.zz();
    const scalar shY  = dF.xy() + dF.yx();
    const scalar shZ  = dF.xz() + dF.zx();

    // Wave strengths (already paired)
    const scalar aP = Up*rho*dF.xx() + (lambda/Up)*trF;
    const scalar aS_y = Us*rho*shY;
    const scalar aS_z = Us*rho*shZ;

    // Roe dissipation (lm)
    fluxLmHat.x() += -0.5*aP;
    fluxLmHat.y() += -0.5*aS_y;
    fluxLmHat.z() += -0.5*aS_z;

    // Roe dissipation (F)
    fluxFHat.xx() += -0.5*(aP)/(rho*Up);
    fluxFHat.yx() += -0.5*(aS_y)/(rho*Us);
    fluxFHat.zx() += -0.5*(aS_z)/(rho*Us);

    // Rotate back
    const vector roeLm = RT & fluxLmHat;
    const tensor roeF = (RT & fluxFHat) & R;

    lmFlux = roeLm * magSf;

    const vector lmC = rho*(roeF & N);
    FFlux = (lmC/rho)*Sf;
}


// ************************************************************************* //
