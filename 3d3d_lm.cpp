#include <getfem/getfem_generic_assembly.h> 
#include <getfem/getfem_assembling.h> 
#include <getfem/getfem_import.h>
#include <getfem/getfem_export.h>   
#include <getfem/getfem_regular_meshes.h>
#include <getfem/getfem_mesh.h>
#include <getfem/getfem_derivatives.h>
#include <getfem/getfem_superlu.h>
#include <getfem/getfem_partial_mesh_fem.h>
#include <getfem/getfem_interpolated_fem.h>
#include <gmm/gmm.h>
#include <gmm/gmm_inoutput.h>
#include <gmm/gmm_iter_solvers.h>
// Standard libraries
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <string>

#include <node.hpp>
#include <math.h>
#include <definitions.hpp>
#include <utilities.hpp>
#include "./muparser/include/muParser.h"

using namespace getfem;


int main(){

/////////////////////////////////////////////////////////////////////
//Importing the parameters of the problem
    ;
//std::string MESH_FILE("/home/federica/Documenti/Federica/FEDE_POLIMI/Software/Stefano/MANworks/transport/src/9_model_error/mesh3D/conforming_mesh8_051_R01.msh");
std::string MESH_FILE("/home/federica/Documenti/Software/Stefano/MANworks/3d3d_lm/cube_nocad.msh");
std::string MESH_TYPE("GT_PK(3,1)");    //GT_PK(n,k) = geometric transformation, n= dim(convexes), k=degree
std::string FEM_TYPEE("FEM_PK(3,1)");   //FEM_PK(n,k) = polynomials on convexes of dim n and of degree k
std::string FEM_TYPEI("FEM_PK(3,1)");
std::string FEM_TYPELM("FEM_PK(3,1)");
std::string IM_TYPE("IM_TETRAHEDRON(6)");
    

/*size_type SIGMA = 6;
size_type OMEGA = 7;
size_type GAMMA = 8;
*/
size_type SIGMA = 70;
size_type OMEGA = 7;
size_type GAMMA =80;

wx = 0.25; //inner x half width
wy = 0.25; //inner y half widht
Wx = 0.5; //outer x half width
Wy = 0.5; //outer y half width
H = 1; //height

// I multiplied by sin (pi*z/H) the second term in ue to have hom bc on top and bottom as in Sigma
//exact solution external domain
std::string ue_ex="(z/H)^2*sin(_pi/H*z)*sin(_pi/wx*x)*sin(_pi/wy*y) + cos(_pi*(x-wx)*(x+wx))*cos(_pi*(y-wy)*(y+wy))*sin(_pi/H*z)";
//exact solution internal domain
std::string ui_ex="(z/H)^2*sin(_pi/H*z)*sin(_pi/wx*x)*sin(_pi/wy*y)";
//source term internal domain
std::string G_EXPR="((_pi/wx)^2 + (_pi/wy)^2) * (z/H)^2*sin(_pi/H*z)*sin(_pi/wx*x)*sin(_pi/wy*y) - sin(_pi/wx*x)*sin(_pi/wy*y)*(2/(H^2)*sin(_pi*z/H)+4*_pi/(H^3)*z*cos(_pi/H*z)-(_pi/H)^2*(z/H)^2*sin(_pi/H*z))";
//source term external domain
std::string F_EXPR="((_pi/wx)^2 + (_pi/wy)^2) * (z/H)^2*sin(_pi/H*z)*sin(_pi/wx*x)*sin(_pi/wy*y) - sin(_pi/wx*x)*sin(_pi/wy*y)*(2/(H^2)*sin(_pi*z/H)+4*_pi/(H^3)*z*cos(_pi/H*z)-(_pi/H)^2*(z/H)^2*sin(_pi/H*z)) + cos(_pi*(y-wy)*(y+wy)) * sin(_pi*z/H)*(-2*_pi* sin(_pi*(x-wx)*(x+wx))+ 4*_pi^2*x^2*cos(_pi*(x-wx)*(x+wx))) + cos(_pi*(x-wx)*(x+wx)) * sin(_pi*z/H)* (-2*_pi* sin(_pi*(y-wy)*(y+wy)) + 4*_pi^2*y^2*cos(_pi*(y-wy)*(y+wy))) + (_pi/H)^2*sin(_pi*z/H)*cos(_pi*(x-wx)*(x+wx))*cos(_pi*(y-wy)*(y+wy))";
//ui-ue
std::string V_EXPR= "- cos(_pi*(x-wx)*(x+wx))*cos(_pi*(y-wy)*(y+wy))*sin(_pi/H*z)";

//////////////////////////////////////////////////////////////////////
//Build the meshes for external and internal domains
    
mesh mesht; //mesh Omega and Sigma
std::cout << "Importing the 3D mesh from gmsh...  "   << std::endl;    
std::string st("gmsh:"+MESH_FILE);
import_mesh(st,mesht);


//////////////////////////////////////////////////////////////////////
//Building the boundaries

std::cout << "Building tissue boundary ..." << std::endl;

size_type face=0;
std::vector< node > BCt;
BCt.clear();

bgeot::pgeometric_trans pgt_t = bgeot::geometric_trans_descriptor(MESH_TYPE);
size_type DIMT = pgt_t->dim();
BCt.reserve(2*DIMT);
// Parse BC data
std::string label_in = "DIR  DIR DIR  DIR  DIR  DIR";
std::string value_in = "0.0  0.0  0.0  0.0  0.0  0.0";

std::vector<std::string> labels = split(label_in, ' ');
std::vector<std::string> values = split(value_in, ' ');
GMM_ASSERT1(labels.size()==2*DIMT, "wrong number of BC labels");
GMM_ASSERT1(values.size()==2*DIMT, "wrong number of BC values");
for (unsigned f=0; f<2*DIMT; ++f) {
    BCt.emplace_back(labels[f], std::stof(values[f]), 0, face+f);
    std::cout << "  face " << f << " : " << BCt.back() << std::endl;
} 


// Build mesht regions
size_type xx=0, yy=1, zz=2;

mesh_region border_faces;
outer_faces_of_mesh(mesht, border_faces);

for (mr_visitor i(border_faces); !i.finished(); ++i) {

    assert(i.is_face());
    base_node un = mesht.normal_of_face_of_convex(i.cv(), i.f());
    un /= gmm::vect_norm2(un);

            if (gmm::abs(un[xx] + 1.0) < 1.0E-7) 	// back
        mesht.region(face+0).add(i.cv(), i.f());
    else if (gmm::abs(un[xx] - 1.0) < 1.0E-7) 	// front
        mesht.region(face+1).add(i.cv(), i.f());
    else if (gmm::abs(un[yy] + 1.0) < 1.0E-7) 	// left
        mesht.region(face+2).add(i.cv(), i.f());
    else if (gmm::abs(un[yy] - 1.0) < 1.0E-7) 	// right
        mesht.region(face+3).add(i.cv(), i.f());
    else if (gmm::abs(un[zz] + 1.0) < 1.0E-7) 	// bottom
        mesht.region(face+4).add(i.cv(), i.f());
    else if (gmm::abs(un[zz] - 1.0) < 1.0E-7) 	// top
        mesht.region(face+5).add(i.cv(), i.f());
            
    }


mesh_region gamma_region;
outer_faces_of_mesh(mesht, mesht.region(OMEGA), gamma_region);

size_type counter_gamma=0;

for (mr_visitor i(gamma_region); !i.finished(); ++i) {

    assert(i.is_face());

    if (!(mesht.region(face+0).is_in(i.cv(),i.f()))&&!(mesht.region(face+1).is_in(i.cv(),i.f()))
        &&!(mesht.region(face+2).is_in(i.cv(),i.f()))&&!(mesht.region(face+3).is_in(i.cv(),i.f()))
        &&!(mesht.region(face+4).is_in(i.cv(),i.f()))&&!(mesht.region(face+5).is_in(i.cv(),i.f()))){ 	// back
        mesht.region(GAMMA).add(i.cv(), i.f());
        counter_gamma+=1;}
    } 
  
mesh_region gamma_region2;
outer_faces_of_mesh(mesht, mesht.region(SIGMA), gamma_region2);

size_type counter_gamma_1=0;

for (mr_visitor i(gamma_region2); !i.finished(); ++i) {

    assert(i.is_face());

    if (!(mesht.region(face+0).is_in(i.cv(),i.f()))&&!(mesht.region(face+1).is_in(i.cv(),i.f()))
    &&!(mesht.region(face+2).is_in(i.cv(),i.f()))&&!(mesht.region(face+3).is_in(i.cv(),i.f()))
    &&!(mesht.region(face+4).is_in(i.cv(),i.f()))&&!(mesht.region(face+5).is_in(i.cv(),i.f()))){ // back
    mesht.region(GAMMA+1).add(i.cv(), i.f());	
    counter_gamma_1 +=1;}
    }
    
std::cout<< "*****nb of faces in gamma " << counter_gamma << "  ****nb of faces in gamma+1 " << counter_gamma_1<<std::endl;

    
    
GMM_ASSERT1(mesht.has_region(GAMMA), "Something wrong happened: I couldn't build the region GAMMA");
GMM_ASSERT1(mesht.has_region(GAMMA+1), "Something wrong happened: I couldn't build the region GAMMA+1");



//////////////////////////////////////////////////////////////////////
//Setting the FEM and the integration method

std::cout << "Setting the FEM and the integration method" << std::endl;

mesh_fem mf_ue(mesht), mf_ui(mesht), mf_p(mesht), mf_p1(mesht);

pfem pf_ue=fem_descriptor(FEM_TYPEE);
pfem pf_ui=fem_descriptor(FEM_TYPEI);
pfem pf_p=fem_descriptor(FEM_TYPELM);
pfem pf_p1=fem_descriptor(FEM_TYPELM);

mf_ue.set_finite_element(mesht.region(OMEGA).index(), pf_ue);
mf_ui.set_finite_element(mesht.region(SIGMA).index(), pf_ui);
mf_p.set_finite_element(mesht.region(GAMMA).index(), pf_p);
mf_p1.set_finite_element(mesht.region(GAMMA+1).index(), pf_p);

mesh_im mimt(mesht);
pintegration_method pim_t = int_method_descriptor(IM_TYPE);
mimt.set_integration_method(mesht.convex_index(), pim_t);

std::cout << "p "<<mf_p.nb_dof() << "   p1  "<< mf_p1.nb_dof() << std::endl;

////////////////////////////////////////////////////////////////////////////////////
//Defining the 3d-3d problem matrix

//total nb of dof 
size_type dof_tot = mf_ue.nb_dof() + mf_ui.nb_dof() + mf_p.nb_dof();

//Global matrix
sparse_matrix_type A(dof_tot, dof_tot); 

//Blocks
sparse_matrix_type Aee(mf_ue.nb_dof(), mf_ue.nb_dof()); 
sparse_matrix_type Aii(mf_ui.nb_dof(), mf_ui.nb_dof());
sparse_matrix_type Bep(mf_ue.nb_dof(), mf_p.nb_dof()); 
sparse_matrix_type Bip(mf_ui.nb_dof(), mf_p.nb_dof()); 

getfem::asm_stiffness_matrix_for_homogeneous_laplacian(Aee,mimt,mf_ue);
getfem::asm_stiffness_matrix_for_homogeneous_laplacian(Aii,mimt,mf_ui);
getfem::asm_mass_matrix (Bep, mimt, mf_ue,mf_p,GAMMA);
gmm::scale(Bep, -1.0);
getfem::asm_mass_matrix(Bip, mimt, mf_ui,mf_p,GAMMA+1);
gmm::add(Aee,  
        gmm::sub_matrix(A, 
                gmm::sub_interval(0, mf_ue.nb_dof()), 
                gmm::sub_interval(0, mf_ue.nb_dof()))); 

gmm::add(Aii,  
        gmm::sub_matrix(A, 
                gmm::sub_interval(mf_ue.nb_dof(), mf_ui.nb_dof()), 
                gmm::sub_interval(mf_ue.nb_dof(), mf_ui.nb_dof())));

gmm::add(Bep,  
        gmm::sub_matrix(A, 
                gmm::sub_interval(0, mf_ue.nb_dof()), 
                gmm::sub_interval(mf_ue.nb_dof() + mf_ui.nb_dof(), mf_p.nb_dof())));

gmm::add(Bip,  
        gmm::sub_matrix(A, 
                gmm::sub_interval(mf_ue.nb_dof(), mf_ui.nb_dof()), 
                gmm::sub_interval(mf_ue.nb_dof() + mf_ui.nb_dof(), mf_p.nb_dof())));

gmm::add(gmm::transposed(Bep),  
        gmm::sub_matrix(A, 
                gmm::sub_interval(mf_ue.nb_dof() + mf_ui.nb_dof(), mf_p.nb_dof()), 
                gmm::sub_interval(0, mf_ue.nb_dof())));

gmm::add(gmm::transposed(Bip),  
        gmm::sub_matrix(A, 
                gmm::sub_interval(mf_ue.nb_dof() + mf_ui.nb_dof(), mf_p.nb_dof()), 
                gmm::sub_interval(mf_ue.nb_dof(), mf_ui.nb_dof())));


///////////////////////////////////////////////////////////////
//Defining the 3d-3d problem RHS
//Global RHS
vector_type F(dof_tot); F.clear();

//source term outer domain
vector_type Fe(mf_ue.nb_dof());
vector_type f_vec(mf_ue.nb_dof()); 
expr = F_EXPR;
interpolation_function(mf_ue, f_vec, f_function); 
getfem::asm_source_term(Fe, mimt, mf_ue, mf_ue, f_vec);

gmm::add(Fe, gmm::sub_vector(F, gmm::sub_interval(0, mf_ue.nb_dof())));

//source term inner domain
vector_type Fi(mf_ui.nb_dof());
vector_type g_vec(mf_ui.nb_dof());
expr=G_EXPR;
interpolation_function(mf_ui, g_vec, f_function);
getfem::asm_source_term(Fi, mimt, mf_ui, mf_ui, g_vec);

gmm::add(Fi, gmm::sub_vector(F, gmm::sub_interval(mf_ue.nb_dof(), mf_ui.nb_dof())));

//source term from ui-ue
vector_type Fp(mf_p.nb_dof());
vector_type v_vec(mf_p.nb_dof());
expr=V_EXPR;
interpolation_function(mf_p, v_vec, f_function);
getfem::asm_source_term(Fp, mimt, mf_p, mf_p, v_vec);

gmm::add(Fp, gmm::sub_vector(F, gmm::sub_interval(mf_ue.nb_dof()+mf_ui.nb_dof(), mf_p.nb_dof())));


//////////////////////////////////////////////////////////////
//Adding boundary conditions
vector_type Dirichlet(mf_ue.nb_dof());
expr=ue_ex;
interpolation_function(mf_ue, Dirichlet, f_function);
for (size_type bc=0; bc < BCt.size(); ++bc)
getfem::assembling_Dirichlet_condition(A, F, mf_ue, BCt[bc].rg, Dirichlet);


//////////////////////////////////////////////////////////////
// Solving the system

    
scalar_type cond;
vector_type U(dof_tot);

size_type restart = 50;
gmm::iteration iter(1e-10);  // iteration object with the max residu
iter.set_noisy(1);               // output of iterations (2: sub-iteration)
iter.set_maxiter(1000);
gmm::identity_matrix PM;
gmm::gmres(A, U, F,PM, restart, iter);

// for(size_type i=0; i<dof_tot; i++)
//     std::cout<< U[i] << std::endl;

gmm::SuperLU_solve(A, U, F, cond);
std::cout << cond << std::endl;

///////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//Test to verify that Bep is assembled well
std::string f=("x^2-2*y^2 -z*y - x*y");
std::string g=("z^2-3*y^2 -4*z*y - x*y");

mesh_fem mf_u(mesht);

pfem pf_u=fem_descriptor("FEM_PK(3,2)");

mf_u.set_finite_element(mesht.region(GAMMA).index(), pf_u);

vector_type fv(mf_u.nb_dof());
interpolation_function(mf_u, fv, f_function);
vector_type gv(mf_u.nb_dof());
interpolation_function(mf_u, gv, f_function);

sparse_matrix_type M(mf_u.nb_dof(), mf_u.nb_dof());
asm_mass_matrix(M, mimt, mf_u, mf_u, GAMMA);
gmm::mult(M, fv, fv); 
scalar_type resu=gmm::vect_sp(f_vec, g_vec);
std::cout << gmm::vect_sp(f_vec, g_vec)<<std::endl;

getfem::ga_workspace workspace;
getfem::size_type nbdofu = mf_u.nb_dof();
getfem::base_vector UU(nbdofu);
workspace.add_fem_variable("u", mf_u, gmm::sub_interval(0, nbdofu), UU);
workspace.add_expression("(X(1)*X(1)-2*X(2)*X(2) -X(3)*X(2) - X(1)*X(2) )* (X(3)*X(3)-3*X(2)*X(2) -4*X(3)*X(2) - X(1)*X(2))", mimt, GAMMA);
// getfem::base_vector L(nbdofu);
// workspace.set_assembled_vector(L);
workspace.assembly(0);

std::cout << workspace.assembled_potential() << std::endl;

/*
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
//This is just a test solving -\Delta u =f in Sigma. It works because we properly enforce the BC on 
//\partial Sigma

sparse_matrix_type M(mf_ui.nb_dof(), mf_ui.nb_dof());
vector_type b(mf_ui.nb_dof());

asm_stiffness_matrix_for_homogeneous_laplacian(M, mimt, mf_ui);


vector_type vec(mf_ui.nb_dof());
expr=G_EXPR;
interpolation_function(mf_ui, vec, f_function);
asm_source_term(b, mimt, mf_ui, mf_ui, vec);

vector_type bc_vector(mf_ui.nb_dof());
expr=ui_ex;
interpolation_function(mf_ui, bc_vector, f_function);

mesh_region sigma_bndary;
outer_faces_of_mesh(mesht, mesht.region(SIGMA), sigma_bndary);
mesht.region(100)=sigma_bndary;
assembling_Dirichlet_condition(M, b, mf_ui, 100, bc_vector );

std::cout <<"solving dummy..."<<std::endl;
gmm::SuperLU_solve(M, UU, b, cond);


mesh_fem mf_uex(mesht);
pfem pf_uex=fem_descriptor("FEM_PK(3,3)");
mf_uex.set_finite_element(mesht.region(SIGMA).index(), pf_uex);
expr=ui_ex;
vector_type vec_uex(mf_uex.nb_dof());
interpolation_function(mf_uex, vec_uex, f_function);


std::vector<scalar_type> L2v(1);
getfem::generic_assembly
	assemL2("Ch=data$1(#1);" "Cex=data$2(#2);"
		"t1=comp(Base(#1).Base(#1))(i,j).Ch(i).Ch(j);"
		"t2=comp(Base(#2).Base(#2))(i,j).Cex(i).Cex(j);"
		"t3=comp(Base(#1).Base(#2))(i,j).Ch(i).Cex(j);"
	"V$1() += t1 + t2 - t3-t3;");
	assemL2.push_mi(mimt);
    assemL2.push_mf(mf_ui);
	assemL2.push_mf(mf_uex);
	assemL2.push_data(UU);
	assemL2.push_data(vec_uex);
	assemL2.push_vec(L2v);
	assemL2.assembly(SIGMA);   //region 0 is the face x=0			

    scalar_type L2_norm=L2v[0];

	L2_norm=sqrt(L2_norm);
	
	cout << "L2 norm error = " << L2_norm << endl;
 */
  
} //end of main


