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

Class
    Foam::localTimeStep

Author
    Oliver Borm
    Aleksandar Jemcov

\*---------------------------------------------------------------------------*/

#include "localTimeStep.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::localTimeStep::localTimeStep
(
    const volVectorField& U,
    const basicThermo& thermo,
    compressible::myTurbulenceModel& turbModel
)
:
    mesh_(U.mesh()),
    U_(U),
    thermo_(thermo),
    turbModel_(turbModel),
    CoDeltaT_
    (
        IOobject
        (
            "CoDeltaT",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("CoDeltaT", dimTime, 0.1)
    ),
    deltaTInvis_
    (
        IOobject
        (
            "deltaTInvis",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("CoDeltaT", dimTime, 0.1)
    )
{


}


void Foam::localTimeStep::update
(
    const scalar maxCo,
    const bool localDt
)
{
    const unallocLabelList& owner = mesh().owner();
    const unallocLabelList& neighbour = mesh().neighbour();

    // Compute characteristic length for each cell
    // Calculated from min face delta coefficient.  HJ, 6/Sep/2012
    // volScalarField deltaX
    // (
    //     IOobject
    //     (
    //         "deltaX",
    //         mesh().time().timeName(),
    //         mesh(),
    //         IOobject::NO_READ,
    //         IOobject::NO_WRITE
    //     ),
    //     mesh(),
    //     dimensionedScalar("great", dimLength, GREAT)
    // );

    // const surfaceScalarField deltaFace
    // (
    //     1/mesh().surfaceInterpolation::deltaCoeffs()
    // );

    // // Compute maximum face area for each cell from the internal faces
    // forAll (owner, facei)
    // {
    //     deltaX[owner[facei]] =
    //         Foam::min(deltaX[owner[facei]], deltaFace[facei]);

    //     deltaX[neighbour[facei]] =
    //         Foam::min(deltaX[neighbour[facei]], deltaFace[facei]);
    // }

    // // Compute maximum face area for each cell from the boundary faces
    // forAll (deltaX.boundaryField(), patchi)
    // {
    //     const fvsPatchScalarField& pDeltaFace =
    //         deltaFace.boundaryField()[patchi];

    //     const fvPatch& p = pDeltaFace.patch();

    //     if (p.coupled())
    //     {
    //         const unallocLabelList& faceCells = p.patch().faceCells();

    //         forAll (pDeltaFace, patchFacei)
    //         {
    //             deltaX[faceCells[patchFacei]] = Foam::min
    //             (
    //                 deltaX[faceCells[patchFacei]],
    //                 pDeltaFace[patchFacei]
    //             );
    //         }
    //     }
    // }

    // store kappa as field
    // need a copy here, because the return value is a <tmp>
    // const volScalarField kappa = thermo_.Cp()/thermo_.Cv();
    // const volScalarField rho = thermo_.rho();
    // const volScalarField p = thermo_.p();

    // // compute the maximum inviscid deltaT
    const volScalarField kappa = thermo_.Cp()/thermo_.Cv();
    volScalarField deltaTInvis0("deltaTInvis0", deltaTInvis_);

    surfaceScalarField amaxSf("amaxSf", 
    mag(fvc::interpolate(U_) & mesh().Sf()) +
    mesh().magSf() * fvc::interpolate(sqrt(thermo_.Cp()/thermo_.Cv()/thermo_.psi())));
    
    scalarField sumAmaxSf(fvc::surfaceSum(amaxSf)().internalField());

   deltaTInvis_.internalField() =  (2*maxCo*mesh().V()) /sumAmaxSf;


        // Update the boundary values
    deltaTInvis_.correctBoundaryConditions();

    // smoothing dt
    scalar localDtSmoothingCoeff
        (
            mesh_.time().controlDict().lookupOrDefault<scalar>
            (
                "localDtSmoothingCoeff",
                0.0
            )
        );
    
        // if (localDtSmoothingCoeff > SMALL)
        // {
        //     deltaTInvis_ +=
        //         localDtSmoothingCoeff
        //     *fvc::laplacian
        //         (
        //             dimensionedScalar("one", dimless, 1.0),
        //             deltaTInvis_
        //         );

        //     deltaTInvis_.correctBoundaryConditions();
        // }
    
    // Info<< "Smoothed flow time scale min/max = "
    //     << gMin(deltaTInvis_.internalField())
    //     << ", " << gMax(deltaTInvis_.internalField()) << endl;

    scalar localDtDampingCoeff
        (
            mesh_.time().controlDict().lookupOrDefault<scalar>
            (
                "localDtDampingCoeff",
                1.0
            )
        );
    	  // Limit rate of change of time scale
     // - reduce as much as required
     // - only increase at a fraction of old time scale
     if
     (
         localDtDampingCoeff < 1.0
      && mesh_.time().timeIndex() > mesh_.time().startTimeIndex() + 1
     )
     {
         deltaTInvis_ =
             deltaTInvis0
            *max(deltaTInvis_/deltaTInvis0, scalar(1) - localDtDampingCoeff);

        //  Info<< "Damped flow time scale min/max = "
        //      << gMin(1/deltaTInvis_.internalField())
        //      << ", " << gMax(1/deltaTInvis_.internalField()) << endl;
     }


    // if (max(turbModel_.muEff()).value() > SMALL)
    // {
    //     // Compute the maximum viscous deltaT
    //     volScalarField deltaTVis
    //     (
    //         ((2*maxCo)*Foam::sqr(deltaX)*thermo_.rho())/turbModel_.muEff()
    //     );
            
    //     //! if the phisicial time detaT() is not fixied this could fails
    //     CoDeltaT_ = min(Foam::min(deltaTVis, deltaTInvis_), (2.0/3.0)* mesh().time().deltaT());
    //     // CoDeltaT_ = Foam::min(deltaTVis, deltaTInvis_);
    // }
    // else
    // {
        CoDeltaT_ = min(deltaTInvis_, (2.0/3.0)* mesh().time().deltaT());
        // CoDeltaT_ = deltaTInvis_;

    // }

    // Info<< "Flow time scale min/max = "
    //     << gMin(CoDeltaT_.internalField())
    //     << ", " << gMax(CoDeltaT_.internalField()) << endl;
    
    if (!localDt)
    {
        CoDeltaT_ = Foam::min(CoDeltaT_);
    }


}


// ************************************************************************* //
