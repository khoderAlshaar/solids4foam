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

#include "makeBasicNumericFlux.H"

// #include "hllcSGFlux.H"
// #include "hllcLMSGFlux.H"
// #include "SGL2RoeFlux.H"
// #include "SGThornberRoeFlux.H"


// #include "roeFlux.H"
// #include "L2RoeFlux.H"
// #include "hlleTNPFlux.H"
// #include "hlleTNPRotatedFlux.H"
// #include "hllcTNPRotatedFlux.H"
// #include "hllcRotatedFlux.H"
// #include "hllcALEFlux.H"
#include "roeALEFlux.H"
#include "L2RoeALEFlux.H"


 
#include "firstOrderLimiter.H" 
#include "BarthJespersenLimiter.H"
#include "VenkatakrishnanLimiter.H" 
#include "Venkatakrishnan5Limiter.H" 
// #include "MinmodLimiter.H"
// #include "SuperbeeLimiter.H"
// #include "VanAlbadaLimiter.H"
#include "WangLimiter.H"
#include "MichalakGoochLimiter.H"
#include "VanLeerLimiter.H"
#include "VanAlbadaLimiter.H"
#include "VenkatakrishnanStableLimiter.H"
// #include "UnlimitedLimiter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * Private Static Data * * * * * * * * * * * * * */

#define makeBasicNumericFluxForAllLimiters(Flux)                              \
makeBasicNumericFlux(Flux, firstOrderLimiter);                                \
makeBasicNumericFlux(Flux, BarthJespersenLimiter);                            \
makeBasicNumericFlux(Flux, VenkatakrishnanLimiter);                       \
makeBasicNumericFlux(Flux, Venkatakrishnan5Limiter);                       \
makeBasicNumericFlux(Flux, WangLimiter);                                    \
makeBasicNumericFlux(Flux, MichalakGoochLimiter);                             \
makeBasicNumericFlux(Flux, VanLeerLimiter);                             \
makeBasicNumericFlux(Flux, VanAlbadaLimiter);                             \
// makeBasicNumericFlux(Flux, VenkatakrishnanStableLimiter);                             \


// makeBasicNumericFluxSG(Flux, SuperbeeLimiter);                                  \
makeBasicNumericFluxSG(Flux, VanAlbadaLimiter);                                 \
// makeBasicNumericFluxSG(Flux, VanLeerLimiter);                                 \
// makeBasicNumericFluxSG(Flux, MinmodLimiter);                                    \
// makeBasicNumericFluxSG(Flux, WangLimiter);                                    \
// makeBasicNumericFluxSG(Flux, MichalakGoochLimiter);                             \
// makeBasicNumericFluxSG(Flux, UnlimitedLimiter);                             \
// makeBasicNumericFlux(Flux, BarthJespersenNewLimiter);                           \

 
// makeBasicNumericFluxForAllLimiters(hllcSGFlux);
// makeBasicNumericFluxForAllLimiters(hllcSGLMFlux);
// makeBasicNumericFluxForAllLimiters(hllcSGLMALEFlux);
// makeBasicNumericFluxForAllLimiters(hllcSGALEFlux);
// makeBasicNumericFluxForAllLimiters(roeSGFlux);
// makeBasicNumericFluxForAllLimiters(hllcSGLMPFlux);
// makeBasicNumericFluxForAllLimiters(SGThornberRoeFlux);
// makeBasicNumericFluxForAllLimiters(SGL2RoeFlux);
// makeBasicNumericFluxForAllLimiters(SGL2RoeALEFlux);
// makeBasicNumericFluxForAllLimiters(SGL2RoeALE2Flux);
// makeBasicNumericFluxForAllLimiters(SGRoeALEFlux);
// makeBasicNumericFluxForAllLimiters(roeFlux);


// makeBasicNumericFluxForAllLimiters(roeFlux);
// makeBasicNumericFluxForAllLimiters(L2RoeFlux);
// makeBasicNumericFluxForAllLimiters(hlleTNPFlux);
// makeBasicNumericFluxForAllLimiters(hllcTNPRotatedFlux);
// makeBasicNumericFluxForAllLimiters(hllcRotatedFlux);
// makeBasicNumericFluxForAllLimiters(hllcALEFlux);
makeBasicNumericFluxForAllLimiters(roeALEFlux);
makeBasicNumericFluxForAllLimiters(L2RoeALEFlux);
// makeBasicNumericFluxForAllLimiters(hlleTNPRotatedFlux);




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
