/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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

#include "numericalFlux.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(numericalFlux, 0);
    defineRunTimeSelectionTable(numericalFlux, state);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::numericalFlux::numericalFlux(Time& runTime)

{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::numericalFlux::~numericalFlux()
{}


// ************************************************************************* //
Foam::autoPtr<Foam::numericalFlux> Foam::numericalFlux::New
(
    Time& runTime,
    const word& region,
    const dynamicFvMesh& mesh,
    volVectorField& lm,
    pointVectorField& lmN,
    volTensorField& F,
    volTensorField& P,
    solidMaterialModel& model,
    operations& op,
    mechanics& mech,
    gradientSchemes& grad
)
{
    // IOdictionary dict
    // (
    //     IOobject
    //     (
    //         "fvSchemes",
    //         mesh.time().caseSystem(),
    //         runTime,
    //         IOobject::MUST_READ,
    //         IOobject::NO_WRITE,
    //         false  // Do not register
    //     )
    // );

    IOdictionary dict
    (
        // If region == "region0" then read from the main case
        // Otherwise, read from the region/sub-mesh directory e.g.
        // constant/fluid or constant/solid
        bool(region == dynamicFvMesh::defaultRegion)
      ? IOobject
        (
            "fvSchemes",
            runTime.caseSystem(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
      : IOobject
        (
            "fvSchemes",
            runTime.caseSystem(),
            region, // using 'local' property of IOobject
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    const dictionary& subDict = dict.subDict("divSchemes").subDict("numericalFlux");

    word name = word(subDict.lookup("flux"));

    // word name = "contact";

    Info<< "Selecting numericFlux " << name << endl;

    stateConstructorTable::iterator cstrIter =
        stateConstructorTablePtr_->find(name);

    if (cstrIter == stateConstructorTablePtr_->end())
    {
        FatalErrorIn("numericalFlux::New(const fvMesh&)")
            << "Unknown numericalFlux type " << name << nl << nl
            << "Valid numericalFlux types are:" << nl
            << stateConstructorTablePtr_->sortedToc() << nl
            << exit(FatalError);
    }

    return autoPtr<numericalFlux>(cstrIter()(runTime, region, mesh, lm, lmN, F, P, model, op, mech, grad));
}