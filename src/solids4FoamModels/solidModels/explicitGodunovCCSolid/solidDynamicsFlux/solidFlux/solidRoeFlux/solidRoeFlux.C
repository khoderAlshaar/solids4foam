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
    const vector& lm_P,
    const vector& lm_M,
    const tensor& F_P,
    const tensor& F_M,
    const tensor& P_P,
    const tensor& P_M,
    // const vector& t_P,
    // const vector& t_M,
    const tensor& R,
    const tensor& RTranspos,
    // const tensor& S_lm,
    // const tensor& S_t,
    const scalar& rho,
    const scalar& lambda,
    const scalar& Up,
    const scalar& Us//,
    // const scalar& magSf,
    // const vector& Sf,
    // const vector& N
) const
{


    
    const vector lm_M_hat = R & lm_M;
    const vector lm_P_hat = R & lm_P;

    
    const tensor F_M_hat = (R & F_M) & RTranspos;
    const tensor F_P_hat = (R & F_P) & RTranspos;

    const tensor P_hat_M = (R & P_M) & RTranspos;
    const tensor P_hat_P = (R & P_P) & RTranspos;

     //difference fields
    const vector    deltaF_lm  = lm_P_hat - lm_M_hat;
    const tensor    deltaF_F   = F_P_hat -  F_M_hat;


    // Central  Flux
    vector flux_lm_P_hat(vector::zero);
    vector flux_lm_M_hat(vector::zero);

    flux_lm_P_hat.x() = P_hat_P.xx();
    flux_lm_P_hat.y() = P_hat_P.yx();
    flux_lm_P_hat.z() = P_hat_P.zx();
    
    flux_lm_M_hat.x() = P_hat_M.xx();
    flux_lm_M_hat.y() = P_hat_M.yx();
    flux_lm_M_hat.z() = P_hat_M.zx();

    const scalar rRho = 1.0 / rho;


    tensor flux_F_M_hat(tensor::zero);
    tensor flux_F_P_hat(tensor::zero);
    
    flux_F_M_hat.xx() = rRho * lm_M_hat.x();
    flux_F_M_hat.yx() = rRho * lm_M_hat.y();
    flux_F_M_hat.zx() = rRho * lm_M_hat.z();
    
    flux_F_P_hat.xx() = rRho * lm_P_hat.x();
    flux_F_P_hat.yx() = rRho * lm_P_hat.y();
    flux_F_P_hat.zx() = rRho * lm_P_hat.z();



    const vector l1_lm(1,0,0); // X-component
    const vector l2_lm(1,0,0); // X-component
    const vector l3_lm(0,1,0); // Y-component
    const vector l4_lm(0,0,1); // Z-component
    const vector l5_lm(0,1,0); // Y-component
    const vector l6_lm(0,0,1); // Z-component

    const scalar rRhoUp = 1.0 / (rho* Up);
    const scalar rRhoUs = 1.0 / (rho* Us);
     
    tensor l1_F(tensor::zero);
    tensor l2_F(tensor::zero);
    tensor l3_F(tensor::zero);
    tensor l4_F(tensor::zero);
    tensor l5_F(tensor::zero);
    tensor l6_F(tensor::zero);

    l1_F.xx() = - rRhoUp;
    l2_F.xx() =   rRhoUp;

    l3_F.yx() = - rRhoUs; // zx-component
    l4_F.zx() = - rRhoUs; // yx-component
    l5_F.yx() =   rRhoUs; // zx-component
    l6_F.zx() =   rRhoUs; // yx-component




    const scalar alpha1 = 0.5*( - Up*rho*deltaF_F.xx() - (lambda/Up)*(deltaF_F.yy()+ deltaF_F.zz() )  +  deltaF_lm.x() );
    const scalar alpha2 = 0.5*(   Up*rho*deltaF_F.xx() + (lambda/Up)*(deltaF_F.yy()+ deltaF_F.zz() )  +  deltaF_lm.x() );

    const scalar alpha3 = 0.5*( - Us* rho * (deltaF_F.xy() + deltaF_F.yx() ) + deltaF_lm.y() );
    const scalar alpha4 = 0.5*( - Us* rho * (deltaF_F.xz() + deltaF_F.zx() ) + deltaF_lm.z() );
    
    const scalar alpha5 = 0.5*(   Us* rho * (deltaF_F.xy() + deltaF_F.yx() ) + deltaF_lm.y() );
    const scalar alpha6 = 0.5*(   Us* rho * (deltaF_F.xz() + deltaF_F.zx() ) + deltaF_lm.z() );

    const vector  sum_lm  = (Up * alpha1 * l1_lm )
                          + (Up * alpha2 * l2_lm )
                          + (Us * alpha3 * l3_lm )
                          + (Us * alpha4 * l4_lm )
                          + (Us * alpha5 * l5_lm )
                          + (Us * alpha6 * l6_lm );

    const tensor sum_F  = (Up * alpha1 * l1_F )
                         + (Up * alpha2 * l2_F )
                         + (Us * alpha3 * l3_F )
                         + (Us * alpha4 * l4_F )
                         + (Us * alpha5 * l5_F )
                         + (Us * alpha6 * l6_F ) ;



                        

    const vector roeFlux_lm_hat = 0.5*  (flux_lm_P_hat + flux_lm_M_hat  + sum_lm);
    const tensor roeFlux_F_hat  = 0.5*  (flux_F_P_hat  + flux_F_M_hat   + sum_F);

    //rotate the system back
    lmFlux = (RTranspos & roeFlux_lm_hat);
    // surfaceTensorField flux_F = RInvers & roeFlux_F_hat & RInversTranspos;
    FFlux = (RTranspos & roeFlux_F_hat) & R;


}

// ************************************************************************* //
