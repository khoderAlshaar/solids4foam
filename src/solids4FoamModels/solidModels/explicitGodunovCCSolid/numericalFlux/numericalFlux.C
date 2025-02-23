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
    const dynamicFvMesh& mesh,
    const volVectorField& lm,
    const volTensorField& F,
    const volTensorField& P,
    solidMaterialModel& model,
    operations& op,
    mechanics& mech
)
{
    IOdictionary dict
    (
        IOobject
        (
            "fvSchemes",
            runTime.system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false  // Do not register
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

    return autoPtr<numericalFlux>(cstrIter()(runTime, mesh, lm, F, P, model, op, mech));
}