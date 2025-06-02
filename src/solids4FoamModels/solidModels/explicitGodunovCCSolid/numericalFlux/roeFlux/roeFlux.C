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

#include "roeFlux.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace numericalFluxs
{
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

    defineTypeNameAndDebug(roeFlux, 0);
    addToRunTimeSelectionTable(numericalFlux, roeFlux, state);






// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

roeFlux::roeFlux
(
    Time& runTime,
    const word& region,
    const dynamicFvMesh& mesh_,
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
    mesh_(mesh_),
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
    N_((Sf_ / mesh_.magSf())),

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
    t_P_((P_P_ & N_)),

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
    interpolate_(mesh_),

    // Cell-averaged linear momentum
    lmR_(interpolate_.surfaceToVol(lmC_)),

    // Local gradient of cell-averaged linear momentum
    lmRgrad_(grad_.localGradient(lmR_, lmC_)),


    lmFlux_
    (
        IOobject
        (
            "lmFlux",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (linearInterpolate(P_) & mesh_.Sf())
    ),
    FFlux_
    (
        IOobject
        (
            "FFlux",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
       ((1/model_.density())*linearInterpolate(lm) * mesh_.Sf())
    ),

    Fx_(op.decomposeTensorX(F_)),
    Fy_(op.decomposeTensorY(F_)),
    Fz_(op.decomposeTensorZ(F_)),

        // compute F gradient
    FxGrad_(grad.gradient(Fx_)),
    FyGrad_(grad.gradient(Fy_)),
    FzGrad_(grad.gradient(Fz_)),
    phi_F_
    (
        IOobject("phi_F", mesh_),
        mesh_,
        dimensionedTensor("phi_F", dimensionSet(0,0,0,0,0,0,0), tensor::zero)
    ),
    F_M_
    (
        IOobject("F_M", mesh_),
        mesh_,
        tensor::I
    ),        
    F_P_ (F_M_),

    T1_
    (
            IOobject("T1_", mesh_),
            mesh_,
            dimensionedVector("T1_", dimensionSet(0,0,0,0,0,0,0), vector::zero)
    ),
    T2_
    (
        IOobject("T2_", mesh_),
        mesh_,
        dimensionedVector("T2_", dimensionSet(0,0,0,0,0,0,0), vector::zero)
    ),
    R_
    (
        IOobject("R", mesh_),
        mesh_,
        tensor::I
    ),
    RTranspos_
    (
        IOobject("RTranspos", mesh_),
        mesh_,
        tensor::I
    ),

    lm_M_hat_(R_ & lm_M_),
    lm_P_hat_(R_ & lm_P_),
    F_M_hat_(R_ & F_M_ & RTranspos_),
    F_P_hat_(R_ & F_P_ & RTranspos_),
    P_hat_P_(R_ & P_P_ & RTranspos_),
    P_hat_M_(R_ & P_M_ & RTranspos_),
    deltaF_lm_(lm_P_hat_ - lm_M_hat_),
    deltaF_F_(F_P_hat_ - F_M_hat_),
    Up_
    (
        IOobject("Up_", mesh_),
        mesh_,
        model.Up()
    ),

    Us_
    (
        IOobject("Us_", mesh_),
        mesh_,
        model.Us()
    ),
    lambda_
    (
        IOobject("lambda_", mesh_),
        mesh_,
        model.lambda()
    ),
    
// Initialize flux fields
    flux_lm_M_hat_
    (
        IOobject("flux_lm_M_hat", mesh_),
        mesh_,
        dimensionedVector("flux_lm_M_hat", dimensionSet(1,-1,-2,0,0,0,0), vector::zero)
    ),

    flux_lm_P_hat_ 
    (
        IOobject("flux_lm_P_hat", mesh_),
        mesh_,
        dimensionedVector("flux_lm_P_hat", dimensionSet(1,-1,-2,0,0,0,0), vector::zero)
    ),

    flux_F_M_hat_
    (
        IOobject("flux_F_M_hat", mesh_), mesh_,
        dimensionedTensor("flux_F_M_hat", dimensionSet(0,1,-1,0,0,0,0), tensor::zero)
    ),

    flux_F_P_hat_
    (
        IOobject("flux_F_P_hat", mesh_),
        mesh_,
        dimensionedTensor("flux_F_P_hat", dimensionSet(0,1,-1,0,0,0,0), tensor::zero)
    ),
        // Initialize eigenvector fields
    l1_lm_
    (
        IOobject("l1_lm", mesh_),
        mesh_,
        vector::zero
    ),
    l2_lm_(l1_lm_),
    l3_lm_(l1_lm_),
    l4_lm_(l1_lm_),
    l5_lm_(l1_lm_),
    l6_lm_(l1_lm_),
    
    l1_F_
    (
        IOobject("l1_F", mesh_),
        mesh_,
        dimensionedTensor("l1_F", dimensionSet(-1,2,1,0,0,0,0), tensor::zero)
    ),

    l2_F_(l1_F_),
    l3_F_(l1_F_),
    l4_F_(l1_F_),
    l5_F_(l1_F_),
    l6_F_(l1_F_),
    
    // Initialize wave strengths
    alpha1_
    (
        IOobject("alpha1", mesh_),
        mesh_,
        dimensionedScalar("alpha1", dimensionSet(1,-2,-1,0,0,0,0), 0.0)
    ),
    alpha2_(alpha1_),
    alpha3_(alpha1_),
    alpha4_(alpha1_),
    alpha5_(alpha1_),
    alpha6_(alpha1_),
    sum_lm_
    (
        IOobject("sum_lm", mesh_),
        mesh_,
        dimensionedVector("flux_lm_M_hat", dimensionSet(1,-1,-2,0,0,0,0), vector::zero)
    ),
    sum_F_
    (
        IOobject("sum_F", mesh_),
        mesh_,
        dimensionedTensor("flux_F_P_hat", dimensionSet(0,1,-1,0,0,0,0), tensor::zero)
    ),
    roeFlux_lm_
    (
        IOobject("roeFlux_lm_", mesh_),
        mesh_,
        dimensionedVector("roeFlux_lm_", dimensionSet(1,-1,-2,0,0,0,0), vector::zero)
    ),
    roeFlux_F_
    (
        IOobject("roeFlux_F_", mesh_),
        mesh_,
        dimensionedTensor("roeFlux_F_", dimensionSet(0,1,-1,0,0,0,0), tensor::zero)
    ),
    roeFlux_lm_hat_
    (
        IOobject("roeFlux_lm_hat", mesh_),
        mesh_,
        dimensionedVector("flux_lm_M_hat", dimensionSet(1,-1,-2,0,0,0,0), vector::zero)
    ),
    roeFlux_F_hat_
    (
        IOobject("roeFlux_F_hat", mesh_),
        mesh_,
        dimensionedTensor("flux_F_M_hat", dimensionSet(0,1,-1,0,0,0,0), tensor::zero)
    )



{

    Info << "Hello from roeFlux constructor" << endl;
   setRotationalMatrix();


//------assign the flux value according to the rotational inveriant----
//-------------------------------------------------------------------------------------
    forAll(flux_lm_P_hat_, facei)
    {
        flux_lm_P_hat_[facei].x() = P_hat_P_[facei].xx();
        flux_lm_P_hat_[facei].y() = P_hat_P_[facei].yx();
        flux_lm_P_hat_[facei].z() = P_hat_P_[facei].zx();
        
        flux_lm_M_hat_[facei].x() = P_hat_M_[facei].xx();
        flux_lm_M_hat_[facei].y() = P_hat_M_[facei].yx();
        flux_lm_M_hat_[facei].z() = P_hat_M_[facei].zx();
        //deformation gradient tensor flux
        flux_F_M_hat_[facei].xx() = (1/rho_.value()) * lm_M_hat_[facei].x();
        flux_F_M_hat_[facei].yx() = (1/rho_.value()) * lm_M_hat_[facei].y();
        flux_F_M_hat_[facei].zx() = (1/rho_.value()) * lm_M_hat_[facei].z();
       
        flux_F_P_hat_[facei].xx() = (1/rho_.value()) * lm_P_hat_[facei].x();
        flux_F_P_hat_[facei].yx() = (1/rho_.value()) * lm_P_hat_[facei].y();
        flux_F_P_hat_[facei].zx() = (1/rho_.value()) * lm_P_hat_[facei].z();
    }


// //!these are constants all the time
forAll(l1_lm_, faceI)
{
    l1_lm_[faceI].x() = 1.0; // x-component

    l2_lm_[faceI].x() = 1.0; // x-component

    l3_lm_[faceI].y() = 1.0; // Y-component

    l4_lm_[faceI].z() = 1.0; // Z-component

    l5_lm_[faceI].y() = 1.0; // Y-component

    l6_lm_[faceI].z() = 1.0; // Z-component


    //assign values to 
    l1_F_[faceI].xx() = -1.0 / (rho_.value()* Up_[faceI]);

    l2_F_[faceI].xx() =  1.0 /  (rho_.value() * Up_[faceI]); // xx-component

    l3_F_[faceI].yx() = -1.0 /  (rho_.value() * Us_[faceI]); // zx-component
    l4_F_[faceI].zx() = -1.0 /  (rho_.value() * Us_[faceI]); // yx-component
    l5_F_[faceI].yx() =  1.0 /  (rho_.value() * Us_[faceI]); // zx-component
    l6_F_[faceI].zx() =  1.0 /  (rho_.value() * Us_[faceI]); // yx-component

}


forAll(alpha1_,facei)
{
    alpha1_[facei] = 0.5*( - Up_[facei] * rho_.value() * deltaF_F_[facei].xx() - (lambda_[facei]/Up_[facei])*(deltaF_F_[facei].yy()+ deltaF_F_[facei].zz() )  +  deltaF_lm_[facei].x() );
    alpha2_[facei] = 0.5*(   Up_[facei] * rho_.value() * deltaF_F_[facei].xx() + (lambda_[facei]/Up_[facei])*(deltaF_F_[facei].yy()+ deltaF_F_[facei].zz() )  +  deltaF_lm_[facei].x() );

    alpha3_[facei] = 0.5*( - Us_[facei]* rho_.value() * (deltaF_F_[facei].xy() +deltaF_F_[facei].yx() ) + deltaF_lm_[facei].y() );
    alpha4_[facei] = 0.5*( - Us_[facei]* rho_.value() * (deltaF_F_[facei].xz() +deltaF_F_[facei].zx() ) + deltaF_lm_[facei].z() );
    
    alpha5_[facei] = 0.5*(   Us_[facei]* rho_.value() * (deltaF_F_[facei].xy() +deltaF_F_[facei].yx() ) + deltaF_lm_[facei].y() );
    alpha6_[facei] = 0.5*(   Us_[facei]* rho_.value() * (deltaF_F_[facei].xz() +deltaF_F_[facei].zx() ) + deltaF_lm_[facei].z() );
}

sum_lm_ = (Up_ * alpha1_ * l1_lm_ )
       + (Up_ * alpha2_ * l2_lm_ )
       + (Us_ * alpha3_ * l3_lm_ )
       + (Us_ * alpha4_ * l4_lm_ )
       + (Us_ * alpha5_ * l5_lm_ )
       + (Us_ * alpha6_ * l6_lm_ );
sum_F_  = (Up_ * alpha1_ * l1_F_ )
       + (Up_ * alpha2_ * l2_F_ )
       + (Us_ * alpha3_ * l3_F_ )
       + (Us_ * alpha4_ * l4_F_ )
       + (Us_ * alpha5_ * l5_F_ )
       + (Us_ * alpha6_ * l6_F_ ) ;

roeFlux_lm_hat_ = 0.5*  (flux_lm_P_hat_ + flux_lm_M_hat_  + sum_lm_);
roeFlux_F_hat_  = 0.5*  (flux_F_P_hat_  + flux_F_M_hat_   + sum_F_);


//rotate the system back
  roeFlux_lm_= (RTranspos_ & roeFlux_lm_hat_);
// surfaceTensorField flux_F = RInvers & roeFlux_F_hat & RInversTranspos;
  roeFlux_F_ = (RTranspos_ & roeFlux_F_hat_) & R_;

lm_b_.correctBoundaryConditions();
t_b_.correctBoundaryConditions();

forAll(mesh_.boundary(), patchi)
{

        forAll(mesh_.boundary()[patchi], facei)
        {
#ifdef OPENFOAM_NOT_EXTEND        
            roeFlux_lm_.boundaryFieldRef()[patchi][facei] =  
                t_b_.boundaryField()[patchi][facei];

            roeFlux_F_.boundaryFieldRef()[patchi][facei] =  (1/rho_.value())*
                ( lm_b_.boundaryField()[patchi][facei] * N_.boundaryField()[patchi][facei]);
#else
            roeFlux_lm_.boundaryField()[patchi][facei] =  
                t_b_.boundaryField()[patchi][facei];

            roeFlux_F_.boundaryField()[patchi][facei] =  (1/rho_.value())*
                ( lm_b_.boundaryField()[patchi][facei] * N_.boundaryField()[patchi][facei]);

#endif
    
        }
}

lmFlux_ = roeFlux_lm_  *mesh_.magSf();

surfaceScalarField N_norm_squared_ = sqr(mag(N_));
lmC_ = (rho_/N_norm_squared_) * (roeFlux_F_ & N_) ;

FFlux_ = (lmC_/rho_)*mesh_.Sf();


}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
roeFlux::~roeFlux()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void roeFlux::setRotationalMatrix()
{
    
// Loop through all faces in the mesh
forAll(N_, facei) 
{
    // // Initialize an arbitrary vector
    vector arbitraryVector(1, 0, 0); // Default arbitrary vector along x-axis

    // Ensure the arbitrary vector is not parallel to the face normal N_[facei]
    if (mag(arbitraryVector & N_[facei]) > 0.999) // If nearly aligned
    {
        arbitraryVector = vector(0, 1, 0); // Switch to y-axis
        if (mag(arbitraryVector & N_[facei]) > 0.999) // If still nearly aligned
        {
            arbitraryVector = vector(0, 0, 1); // Use z-axis as a last resort
        }
    }

    // Step 3: Project arbitraryVector onto the plane tangent to N_[facei]
    T1_[facei] = arbitraryVector - (arbitraryVector & N_[facei]) * N_[facei];

    // Find an orthogonal vector to N_[facei]
    // T1_[facei] = op.findOrthogonal(N_[facei]);

    // Normalize T1_ to make it a unit vector
    scalar T1_Mag = mag(T1_[facei]);
    if (T1_Mag > SMALL) // Avoid division by zero
    {
        T1_[facei] /= T1_Mag;
    }
    else
    {
        FatalErrorInFunction << "Zero-length T1_ encountered at face " << facei << abort(FatalError);
    }

    // Second tangential vector (T2_): Cross product of N and T1_
    T2_[facei] = N_[facei] ^ T1_[facei];

    // Normalize T2_ to make it a unit vector
    scalar T2_Mag = mag(T2_[facei]);
    if (T2_Mag > SMALL) // Avoid division by zero
    {
        T2_[facei] /= T2_Mag;
    }
    else
    {
        FatalErrorInFunction << "Zero-length T2_ encountered at face " << facei << abort(FatalError);
    }

    // Verify orthogonality of T1_, T2_, and N
    scalar tolerance = 1e-6;

    if (mag(N_[facei] & T1_[facei]) > tolerance) 
    {
        FatalErrorInFunction << "N and T1_ are not orthogonal at face " << facei 
                             << ". Dot product: " << (N_[facei] & T1_[facei]) << abort(FatalError);
    }

    if (mag(N_[facei] & T2_[facei]) > tolerance) 
    {
        FatalErrorInFunction << "N and T2_ are not orthogonal at face " << facei 
                             << ". Dot product: " << (N_[facei] & T2_[facei]) << abort(FatalError);
    }

    if (mag(T1_[facei] & T2_[facei]) > tolerance) 
    {
        FatalErrorInFunction << "T1_ and T2_ are not orthogonal at face " << facei 
                             << ". Dot product: " << (T1_[facei] & T2_[facei]) << abort(FatalError);
    }

    // Now T1_ and T2_ are guaranteed orthogonal to N and each other
}
//these are constant values
    forAll(R_, facei)
    {
        R_[facei].xx() = N_[facei].x();
        R_[facei].xy() = N_[facei].y();
        R_[facei].xz() = N_[facei].z();
        R_[facei].yx() = T1_[facei].x();
        R_[facei].yy() = T1_[facei].y();
        R_[facei].yz() = T1_[facei].z();
        R_[facei].zx() = T2_[facei].x();
        R_[facei].zy() = T2_[facei].y();
        R_[facei].zz() = T2_[facei].z();
    }

RTranspos_ = R_.T();

}



void roeFlux::reconstruction()
{

    // Cell gradients
    lmGrad_ = grad_.gradient(lm_);
    grad_.gradient(P_, PxGrad_, PyGrad_, PzGrad_);
    grad_.gradient(F_, FxGrad_, FyGrad_, FzGrad_);


    // Reconstruction
    grad_.reconstruct(lm_, lmGrad_, lm_M_, lm_P_);
    grad_.reconstruct(P_, PxGrad_, PyGrad_, PzGrad_, P_M_, P_P_);
    grad_.reconstruct(F_, FxGrad_, FyGrad_, FzGrad_, F_M_, F_P_);

    t_M_ = P_M_ & N_;
    t_P_ = P_P_ & N_;

    // Riemann solver
    S_lm_ = mech_.Smatrix_lm();
    S_t_ = mech_.Smatrix_t();
}


void roeFlux::curllFreeAlgorithm()
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

void roeFlux::computeFlux()
{
    // Info << "Compute flux using Contact flux" << endl;

    reconstruction();

// Acoustic Riemann solver
#ifdef OPENFOAM_NOT_EXTEND        
S_lm_.oriented() = false;
S_t_.oriented() = false;
t_M_.oriented() = false;
t_P_.oriented() = false;
lm_P_.oriented() = false;
lm_M_.oriented() = false;

F_P_.oriented() = false;
F_M_.oriented() = false;
 
P_P_.oriented() = false;
P_M_.oriented() = false;
#endif

 lm_M_hat_ = R_ & lm_M_;
 lm_P_hat_ = R_ & lm_P_;



 F_M_hat_ = (R_ & F_M_) & RTranspos_;
 F_P_hat_ = (R_ & F_P_) & RTranspos_;


 P_hat_M_ = (R_ & P_M_) & RTranspos_;
 P_hat_P_ = (R_ & P_P_) & RTranspos_;

 //difference fields
deltaF_lm_  = lm_P_hat_ - lm_M_hat_;
deltaF_F_   = F_P_hat_ -  F_M_hat_;

    forAll(flux_lm_P_hat_, facei)
    {
        flux_lm_P_hat_[facei].x() = P_hat_P_[facei].xx();
        flux_lm_P_hat_[facei].y() = P_hat_P_[facei].yx();
        flux_lm_P_hat_[facei].z() = P_hat_P_[facei].zx();
        
        flux_lm_M_hat_[facei].x() = P_hat_M_[facei].xx();
        flux_lm_M_hat_[facei].y() = P_hat_M_[facei].yx();
        flux_lm_M_hat_[facei].z() = P_hat_M_[facei].zx();
        //deformation gradient tensor flux
        flux_F_M_hat_[facei].xx() = (1/rho_.value()) * lm_M_hat_[facei].x();
        flux_F_M_hat_[facei].yx() = (1/rho_.value()) * lm_M_hat_[facei].y();
        flux_F_M_hat_[facei].zx() = (1/rho_.value()) * lm_M_hat_[facei].z();
       
        flux_F_P_hat_[facei].xx() = (1/rho_.value()) * lm_P_hat_[facei].x();
        flux_F_P_hat_[facei].yx() = (1/rho_.value()) * lm_P_hat_[facei].y();
        flux_F_P_hat_[facei].zx() = (1/rho_.value()) * lm_P_hat_[facei].z();
    }

 //eigenvalues and eigenvectors are constants no needd to update them again for a linear system


forAll(alpha1_,facei)
{
    alpha1_[facei] = 0.5*( - Up_[facei] * rho_.value() * deltaF_F_[facei].xx() - (lambda_[facei]/Up_[facei])*(deltaF_F_[facei].yy()+ deltaF_F_[facei].zz() )  +  deltaF_lm_[facei].x() );
    alpha2_[facei] = 0.5*(   Up_[facei] * rho_.value() * deltaF_F_[facei].xx() + (lambda_[facei]/Up_[facei])*(deltaF_F_[facei].yy()+ deltaF_F_[facei].zz() )  +  deltaF_lm_[facei].x() );

    alpha3_[facei] = 0.5*( - Us_[facei]* rho_.value() * (deltaF_F_[facei].xy() +deltaF_F_[facei].yx() ) + deltaF_lm_[facei].y() );
    alpha4_[facei] = 0.5*( - Us_[facei]* rho_.value() * (deltaF_F_[facei].xz() +deltaF_F_[facei].zx() ) + deltaF_lm_[facei].z() );
    
    alpha5_[facei] = 0.5*(   Us_[facei]* rho_.value() * (deltaF_F_[facei].xy() +deltaF_F_[facei].yx() ) + deltaF_lm_[facei].y() );
    alpha6_[facei] = 0.5*(   Us_[facei]* rho_.value() * (deltaF_F_[facei].xz() +deltaF_F_[facei].zx() ) + deltaF_lm_[facei].z() );
}

sum_lm_ = (Up_ * alpha1_ * l1_lm_ )
       + (Up_ * alpha2_ * l2_lm_ )
       + (Us_ * alpha3_ * l3_lm_ )
       + (Us_ * alpha4_ * l4_lm_ )
       + (Us_ * alpha5_ * l5_lm_ )
       + (Us_ * alpha6_ * l6_lm_ );
sum_F_  = (Up_ * alpha1_ * l1_F_ )
       + (Up_ * alpha2_ * l2_F_ )
       + (Us_ * alpha3_ * l3_F_ )
       + (Us_ * alpha4_ * l4_F_ )
       + (Us_ * alpha5_ * l5_F_ )
       + (Us_ * alpha6_ * l6_F_ ) ;

roeFlux_lm_hat_ = 0.5*  (flux_lm_P_hat_ + flux_lm_M_hat_  + sum_lm_);
roeFlux_F_hat_  = 0.5*  (flux_F_P_hat_  + flux_F_M_hat_   + sum_F_);


//rotate the system back
  roeFlux_lm_= (RTranspos_ & roeFlux_lm_hat_);
// surfaceTensorField flux_F = RInvers & roeFlux_F_hat & RInversTranspos;
  roeFlux_F_ = (RTranspos_ & roeFlux_F_hat_) & R_;

lm_b_.correctBoundaryConditions();
t_b_.correctBoundaryConditions();

forAll(mesh_.boundary(), patchi)
{

        forAll(mesh_.boundary()[patchi], facei)
        {
#ifdef OPENFOAM_NOT_EXTEND 
            roeFlux_lm_.boundaryFieldRef()[patchi][facei] =  
                t_b_.boundaryField()[patchi][facei];

            roeFlux_F_.boundaryFieldRef()[patchi][facei] =  (1/rho_.value())*
                ( lm_b_.boundaryField()[patchi][facei] * N_.boundaryField()[patchi][facei]);
#else
            roeFlux_lm_.boundaryField()[patchi][facei] =  
                t_b_.boundaryField()[patchi][facei];

            roeFlux_F_.boundaryField()[patchi][facei] =  (1/rho_.value())*
                ( lm_b_.boundaryField()[patchi][facei] * N_.boundaryField()[patchi][facei]);
#endif
        }
}

lmFlux_ = roeFlux_lm_  *mesh_.magSf();//correct

surfaceScalarField N_norm_squared_ = sqr(mag(N_));
lmC_ = (rho_/N_norm_squared_) * (roeFlux_F_ & N_) ;

curllFreeAlgorithm();

// roeFlux_F_ = (1/rho_)* (lmC * N);

FFlux_ = (lmC_/rho_)*mesh_.Sf();

    
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace numericalFluxs

} // End namespace Foam

// ************************************************************************* //
