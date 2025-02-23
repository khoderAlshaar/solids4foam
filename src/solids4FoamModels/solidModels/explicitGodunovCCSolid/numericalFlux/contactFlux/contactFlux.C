/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "contactFlux.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace numericalFluxs
{
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

    defineTypeNameAndDebug(contactFlux, 0);
    addToRunTimeSelectionTable(numericalFlux, contactFlux, state);






// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

contactFlux::contactFlux
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
:
    numericalFlux(runTime), //! later on I can change this by moving all these to the numericalFluc class and make access functions
    mesh_(mesh),
    lm_(lm),
    lmN_(lmN),
    F_(F),
    P_(P),
    model_(model),
    op_(op),
    mech_(mech),
    grad_(grad),
//----------------
    magSf_(mesh_.magSf()),
    Sf_(mesh_.Sf()),
    // Creating mesh normal fields
    N_((Sf_ / mesh_.magSf()).ref()),

    rho_(model_.density()),

    lmGrad_(grad_.gradient(lm_)),

    Px_(op_.decomposeTensorX(P_)),
    Py_(op_.decomposeTensorY(P_)),
    Pz_(op_.decomposeTensorZ(P_)),

    PxGrad_(grad_.gradient(Px_)),
    PyGrad_(grad_.gradient(Py_)),
    PzGrad_(grad_.gradient(Pz_)),

    // Reconstruction of linear momentum
    lm_M_(
        IOobject("lm_M", mesh_),
        mesh_,
        dimensionedVector("lm_M", lm_.dimensions(), vector::zero)
    ),

    lm_P_(lm_M_),

    // Reconstruction of PK1 stresses
    P_M_(
        IOobject("P_M", mesh_),
        mesh_,
        dimensionedTensor("P_M", P_.dimensions(), tensor::zero)
    ),
    P_P_(P_M_),

    // Reconstruction of traction
    t_M_
    (
        IOobject("t_M", mesh_),
        P_M_ & N_
    ),
    t_P_((P_P_ & N_).ref()),

    S_lm_(mech_.Smatrix_lm()),
    S_t_(mech_.Smatrix_t()),

    // Contact traction
    tC_(t_M_),
    // Contact linear momentum
    lmC_(lm_M_),
    t_b_
    (
        IOobject
        (
            "t_b",
            runTime.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    lm_b_
    (
        IOobject
        (
            "lm_b",
            runTime.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    phi_lm_
    (
        IOobject("phi_lm", mesh_),
        mesh_,
        dimensionedVector("phi_lm", dimensionSet(0,0,0,0,0,0,0), vector::zero)
    ),
    phi_P_
    (
        IOobject("phi_P", mesh_),
        mesh_,
        dimensionedTensor("phi_P", dimensionSet(0,0,0,0,0,0,0), tensor::zero)
    ),

        // Constrained class
    interpolate_(mesh),

    // Cell-averaged linear momentum
    lmR_(interpolate_.surfaceToVol(lmC_)),

    // Local gradient of cell-averaged linear momentum
    lmRgrad_(grad_.localGradient(lmR_, lmC_)),


    lmFlux_
    (
        IOobject
        (
            "lmFlux",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (linearInterpolate(P_) & mesh.Sf())
    ),
    FFlux_
    (
        IOobject
        (
            "FFlux",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
       ((1/model_.density())*linearInterpolate(lm) * mesh.Sf())
    )

{

    Info << "Hello from contactFlux constructor" << endl;
   
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
contactFlux::~contactFlux()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void contactFlux::reconstruction()
{

    // Cell gradients
    lmGrad_ = grad_.gradient(lm_);
    grad_.gradient(P_, PxGrad_, PyGrad_, PzGrad_);


    // Reconstruction
    grad_.reconstruct(lm_, lmGrad_, lm_M_, lm_P_);
    grad_.reconstruct(P_, PxGrad_, PyGrad_, PzGrad_, P_M_, P_P_);

    t_M_ = P_M_ & N_;
    t_P_ = P_P_ & N_;

    // Riemann solver
    S_lm_ = mech_.Smatrix_lm();
    S_t_ = mech_.Smatrix_t();
}


void contactFlux::curllFreeAlgorithm()
{
    
    lmR_ = interpolate_.surfaceToVol(lmC_);
    lmRgrad_ = grad_.localGradient(lmR_, lmC_);
    interpolate_.volToPoint(lmR_, lmRgrad_, lmN_);
    #include "strongBCs.H"
    lmN_.correctBoundaryConditions();

    // Info << "lmN_"  << lmN_<<endl;

    // Constrained fluxes
    lmC_ = interpolate_.pointToSurface(lmN_);

}

void contactFlux::computeFlux()
{
    // Info << "Compute flux using Contact flux" << endl;

    reconstruction();

// Acoustic Riemann solver
S_lm_.oriented() = false;
S_t_.oriented() = false;
t_M_.oriented() = false;
t_P_.oriented() = false;
lm_P_.oriented() = false;
lm_M_.oriented() = false;

tC_ = 0.5*(t_M_+t_P_) + (0.5*S_lm_ & (lm_P_ - lm_M_));
lmC_ = 0.5*(lm_M_+lm_P_) + (0.5*S_t_ & (t_P_ - t_M_));


// Compute boundary values
lm_b_.correctBoundaryConditions();
t_b_.correctBoundaryConditions();


// if (Pstream::parRun())
// {
//     op_.decomposeTensor(P_, Px_, Py_, Pz_);
//     n_ = mech_.spatialNormal(F_);
// }


forAll(mesh_.boundary(), patchi)
{
    // // Riemann solver for inter-processor boundaries
    // if (mesh_.boundary()[patchi].coupled())
    // {
    //     const vectorField lm_nei(
    //       lm_.boundaryField()[patchi].patchNeighbourField());

    //     const tensorField P_nei(
    //       P_.boundaryField()[patchi].patchNeighbourField());

    //     const vectorField Px_nei(
    //       Px_.boundaryField()[patchi].patchNeighbourField());

    //     const vectorField Py_nei(
    //       Py_.boundaryField()[patchi].patchNeighbourField());

    //     const vectorField Pz_nei(
    //       Pz_.boundaryField()[patchi].patchNeighbourField());

    //     const tensorField lmGrad_nei(
    //       lmGrad_.boundaryField()[patchi].patchNeighbourField());

    //     const tensorField PxGrad_nei(
    //       PxGrad_.boundaryField()[patchi].patchNeighbourField());

    //     const tensorField PyGrad_nei(
    //       PyGrad_.boundaryField()[patchi].patchNeighbourField());

    //     const tensorField PzGrad_nei(
    //       PzGrad_.boundaryField()[patchi].patchNeighbourField());

    //     const vectorField C_nei(
    //       C_.boundaryField()[patchi].patchNeighbourField());

    //     const scalarField Up_nei(
    //       Up_.boundaryField()[patchi].patchNeighbourField());

    //     const scalarField Us_nei(
    //       Us_.boundaryField()[patchi].patchNeighbourField());

    //     forAll(mesh_.boundary()[patchi], facei)
    //     {
    //         const label& bCell =
    //             mesh_.boundaryMesh()[patchi].faceCells()[facei];

    //         const vector& Cf = mesh_.Cf().boundaryField()[patchi][facei];

    //         const vector& lm_M =
    //             lm_[bCell] + (lmGrad_[bCell] & (Cf - C_[bCell]));

    //         const vector& lm_P =
    //             lm_nei[facei] + (lmGrad_nei[facei] & (Cf - C_nei[facei]));

    //         const vector& Px_M =
    //             Px_[bCell] + (PxGrad_[bCell] & (Cf - C_[bCell]));

    //         const vector& Px_P =
    //                 Px_nei[facei] + (PxGrad_nei[facei] & (Cf - C_nei[facei]));

    //         const vector& Py_M =
    //             Py_[bCell] + (PyGrad_[bCell] & (Cf - C_[bCell]));

    //         const vector& Py_P =
    //             Py_nei[facei] + (PyGrad_nei[facei] & (Cf - C_nei[facei]));

    //         const vector& Pz_M =
    //             Pz_[bCell] + (PzGrad_[bCell] & (Cf - C_[bCell]));

    //         const vector& Pz_P =
    //             Pz_nei[facei] + (PzGrad_nei[facei] & (Cf - C_nei[facei]));

    //         const tensor& P_M = tensor(Px_M, Py_M, Pz_M);
    //         const tensor& P_P = tensor(Px_P, Py_P, Pz_P);

    //         const scalar Up = (Up_[bCell] + Up_nei[facei])/2.0;
    //         const scalar Us = (Us_[bCell] + Us_nei[facei])/2.0;

    //         const vector& N = N_.boundaryField()[patchi][facei];
    //         const vector& n = n_.boundaryField()[patchi][facei];
    //         const tensor S_lm = (Up*n*n) + (Us*(I-(n*n)));
    //         const tensor S_t = ((n*n)/Up) + ((I-(n*n))/Us);

    //         tC_.boundaryFieldRef()[patchi][facei] =
    //             0.5*((P_M + P_P) & N) + (0.5*S_lm & (lm_P - lm_M));

    //         lmC_.boundaryFieldRef()[patchi][facei] =
    //             0.5*(lm_M + lm_P) + 0.5*(S_t & ((P_P - P_M) & N));
    //     }
    // }

    // Apply boundary conditions
    // else
    // {
        forAll(mesh_.boundary()[patchi], facei)
        {
            lmC_.boundaryFieldRef()[patchi][facei] =
                lm_b_.boundaryField()[patchi][facei];

            tC_.boundaryFieldRef()[patchi][facei] =
                t_b_.boundaryField()[patchi][facei];
        }
    // }
}

    curllFreeAlgorithm();

    lmFlux_ = tC_*magSf_;
    FFlux_ = (lmC_/rho_)*Sf_ ;

    
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace numericalFluxs

} // End namespace Foam

// ************************************************************************* //
