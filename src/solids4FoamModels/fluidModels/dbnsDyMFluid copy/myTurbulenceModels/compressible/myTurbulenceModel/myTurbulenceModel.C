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

#include "myTurbulenceModel.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(myTurbulenceModel, 0);
defineRunTimeSelectionTable(myTurbulenceModel, myTurbulenceModel);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

myTurbulenceModel::myTurbulenceModel
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const basicThermo& thermophysicalModel,
    const word& myTurbulenceModelName
)
:
    regIOobject
    (
        IOobject
        (
            myTurbulenceModelName,
            U.time().constant(),
            U.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    runTime_(U.time()),
    mesh_(U.mesh()),

    rho_(rho),
    U_(U),
    phi_(phi),
    thermophysicalModel_(thermophysicalModel),

    y_(mesh_)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<myTurbulenceModel> myTurbulenceModel::New
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const basicThermo& thermophysicalModel,
    const word& myTurbulenceModelName
)
{
    word modelName;

    // Enclose the creation of the dictionary to ensure it is deleted
    // before the myTurbulenceModel is created otherwise the dictionary is
    // entered in the database twice
    {
        IOdictionary dict
        (
            IOobject
            (
                "turbulenceProperties",
                U.time().constant(),
                U.db(),
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        );

        dict.lookup("simulationType") >> modelName;
    }

    Info<< "Selecting turbulence model type " << modelName << endl;

    myTurbulenceModelConstructorTable::iterator cstrIter =
        myTurbulenceModelConstructorTablePtr_->find(modelName);

    if (cstrIter == myTurbulenceModelConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "myTurbulenceModel::New(const volScalarField&, "
            "const volVectorField&, const surfaceScalarField&, "
            "basicThermo&)"
        )   << "Unknown myTurbulenceModel type " << modelName
            << endl << endl
            << "Valid myTurbulenceModel types are :" << endl
            << myTurbulenceModelConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<myTurbulenceModel>
    (
        cstrIter()(rho, U, phi, thermophysicalModel, myTurbulenceModelName)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar myTurbulenceModel::yPlusLam(const scalar kappa, const scalar E) const
{
    scalar ypl = 11.0;

    for (int i=0; i<10; i++)
    {
        ypl = log(E*ypl)/kappa;
    }

    return ypl;
}


tmp<volScalarField> myTurbulenceModel::rhoEpsilonEff() const
{
    tmp<volTensorField> tgradU = fvc::grad(U_);
    return mu()*(tgradU() && dev(twoSymm(tgradU()))) + rho_*epsilon();
}


void myTurbulenceModel::correct()
{
    if (mesh_.changing())
    {
        y_.correct();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
