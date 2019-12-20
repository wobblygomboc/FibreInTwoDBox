//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented,
//LIC// multi-physics finite-element library, available
//LIC// at http://www.oomph-lib.org.
//LIC//
//LIC//    Version 1.0; svn revision $LastChangedRevision: 1097 $
//LIC//
//LIC// $LastChangedDate: 2015-12-17 11:53:17 +0000 (Thu, 17 Dec 2015) $
//LIC//
//LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
//LIC//
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC//
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC//
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC//
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC//
//LIC//====================================================================
//Driver for 1D fibre moving in a 2D box -- meshed with triangle

// forces the use of a finite-differenced jacobian, for debug
#define USE_FD_JACOBIAN

//Generic includes
#include "generic.h"

// Navier-Stokes equations and Taylor-Hood elements
#include "navier_stokes.h"

// The mesh
#include "meshes/triangle_mesh.h"

// wrapper for locate zeta to allow us to visualise along a given line
#include "generic/line_visualiser.h"

// classes which contain the machinery for the subtraction of singularities
#include "navier_stokes_sing_face_element.h"

// generic class for a plate as a geometric object, and specific cases,
// a spherical plate and a flat plate
# include "plate_as_geometric_object.h"

// the solution for the flow around a semi-infinite plate at x>0 
#include "moffat_solution.h"

// ********************************************************************************

// define DEBUG for additional output
#define xDEBUG

using namespace std;

using namespace oomph;

namespace Additional_Maths_Functions
{
  double atan2pi(const double y, const double x)
  {
    // Polar angle
    double theta = atan2(y,x);

    // prevent atan2 negative angle fuckery that causes a discontinuity at theta=pi
    if (y < 0.0)
    {
      theta += 2.0 * MathematicalConstants::Pi;
    }

    return theta;
  }

  // cosec, for mathematica output
  double Csc(const double& x)
  {
    return 1.0/sin(x);
  }
  
  // sign function
  template <typename T>
  int sgn(T val)
  {
    return (T(0) < val) - (val < T(0));
  }

  int delta(const int& i, const int& j)
  {
    return (i == j);
  }
}

//==start_of_namespace==============================
/// Namespace for physical parameters
//==================================================
namespace Global_Physical_Variables
{
  // Results output directory
  string dir = "RESLT";
  
  // hacky
  Plate* Plate_pt = 0;

  // fake singular amplitude for debugging
  // (when bypassing proper calculation of the amplitude)
  double singular_amplitude_for_debug = 0;
  
  // Dimensionless domain values
  // ---------------------------
  double L_edge_x      = -0.5;
  double R_edge_x      =  0.5;
  double top_edge_y    =  0.25;
  double bottom_edge_y = -0.25;
  double plate_radius  =  0.25;

  double plate_radius_of_curvature = 100;
  
  // standard element area for the bulk 
  double element_area = 0.5;

  // high resolution element area for the regions near the plate edges
  double hi_res_element_area = 1e-2; // 1e-6;

  // flag to switch the high-res regions on/off
  bool include_hi_res_regions = true;
  
  // double high_res_region_radius = 0.1;  

  // angular coordinate along plate where the circles bounding the high-resolution regions
  // intersect. 0.17rad ~10deg
  double high_res_region_zeta = 0.1;
  
  /// Reynolds number
  double Re = 0.0; 

  // rigid-body velocity of the boundaries (u_x, u_y)
  Vector<double> boundary_velocity(2, 0.0);
  
  // any rotational velocity the plate/boundary has, positive is anticlockwise (i.e. right-hand rule)
  double angular_velocity_plate    = 0.0;
  double angular_velocity_boundary = 0.0;
  
  // the point about which the plate rotates
  Vector<double> centre_of_rotation(2, 0.0);

  // rigid-body velocity of the plate (u_x, u_y)
  Vector<double> plate_velocity(2, 0.0);
  
  // boundaries
  enum
  {
    First_outer_boundary            = 0,
    Outer_boundary_upper            = 0,
    Outer_boundary_lower            = 1,
    Last_outer_boundary,                  // for looping

    First_inner_boundary            = 2,
    Inner_region_boundary_left      = 2,
    Inner_plate_boundary            = 3,
    Inner_region_boundary_right     = 4,
    Inner_region_boundary_right_con = 5,
    Inner_hi_res_region_left_upper  = 6,
    Inner_hi_res_region_left_lower  = 7,
    Inner_hi_res_region_right_upper = 8,
    Inner_hi_res_region_right_connecting = 9,
    Inner_hi_res_region_right_lower = 10,
    Last_inner_boundary
  };

  // this will be assigned once Triangle has assembled the mesh
  // and will change if the high-res regions aren't being included
  unsigned Inner_plate_boundary_upper = 11;
  
  // regions
  enum
  {
    Region_upper              = 1, 
    Region_hi_res_left_upper  = 2,
    Region_hi_res_left_lower  = 3,
    Region_hi_res_right_upper = 4,
    Region_hi_res_right_lower = 5
  };
 
  // Bit of a hack but it facilitates reuse...
  #include "unstructured_moffatt_mesh_3_inner_boundaries.h"

  Vector<double> u_BC(unsigned boundary_id)
  {
    Vector<double> u_bc(2,0.0);
    
    // in this instance, we're just doing no-slip & no-penetration everywhere
    
    return u_bc;
  }
} // end_of_namespace

//==start_of_namespace==============================
/// Namespace for analytic functions
//==================================================
namespace Analytic_Functions
{  
  /// \short Function to convert 2D Polar derivatives (du/dr, du/dtheta, dv/dr, dv/dtheta)
  // to Cartesian derivatives (dux/dx, dux/dy, duy/dx, duy/dy)
  DenseMatrix<double> polar_to_cartesian_derivatives_2d(const DenseMatrix<double> grad_u_polar,
							const Vector<double> u_polar,
							const double r, double theta)
  {
    // shorthand for polar components
    double u = u_polar[0];
    double v = u_polar[1];

    double du_dr     = grad_u_polar(0,0);
    double du_dtheta = grad_u_polar(0,1);
    double dv_dr     = grad_u_polar(1,0);
    double dv_dtheta = grad_u_polar(1,1);
       
    // output cartesian tensor du_i/dx_j
    DenseMatrix<double> du_dx(2, 2, 0.0);

    // dux_dx
    du_dx(0,0) = cos(theta) * (du_dr*cos(theta) - dv_dr*sin(theta))
      -(1.0/r)*sin(theta) * (du_dtheta*cos(theta) - u*sin(theta) -
    			     dv_dtheta*sin(theta) - v*cos(theta));
  
    // dux_dy
    du_dx(0,1) = sin(theta) * (du_dr*cos(theta) - dv_dr*sin(theta))
      +(1.0/r)*cos(theta) * (du_dtheta*cos(theta) - u*sin(theta) -
    			     dv_dtheta*sin(theta) - v*cos(theta));

    // duy_dx
    du_dx(1,0) = cos(theta) * (du_dr*sin(theta) + dv_dr*cos(theta))
      -(1.0/r)*sin(theta) * (du_dtheta*sin(theta) + u*cos(theta) +
    			     dv_dtheta*cos(theta) - v*sin(theta));

    // duy_dy
    du_dx(1,1) = sin(theta) * (du_dr*sin(theta) + dv_dr*cos(theta))
      +(1.0/r)*cos(theta) * (du_dtheta*sin(theta) + u*cos(theta) +
    			     dv_dtheta*cos(theta) - v*sin(theta));

    return du_dx;
  }

  /// \short Newtonian stress tensor
  DenseMatrix<double> newtonian_stress(const DenseMatrix<double>& strain_rate,
				       const double& p)
  {
    // \tau_{ij}
    DenseMatrix<double> stress(2,2);
    
    for(unsigned i=0; i<2; i++)
    {
      for(unsigned j=0; j<2; j++)
      {
	// Newtonian constitutive relation
	stress(i,j) = -p*Additional_Maths_Functions::delta(i,j) + 2.0*strain_rate(i,j);
      }
    }

    return stress;
  }
  
  /// \short "Singular" function and gradient
  void singular_fct_and_gradient(const Vector<double>& x,
				 Vector<double>& u, DenseMatrix<double>& du_dx,
				 int boundary_id=-1)
  {
    // big number to numerically represent infinity; value of 1/r at 1/2 of
    // the element lengthscale
    const double infinity = 2.0 / sqrt(2.0 * Global_Physical_Variables::element_area);
    
    // this coresponds to pure normal flow
    double A = 0;
    double B = -1;

    // location of plate starting edge   
    Vector<double> zeta_start(1);
    zeta_start[0] = Global_Physical_Variables::Plate_pt->zeta_start();

    Vector<double> left_coords(2);
    Global_Physical_Variables::Plate_pt->position(zeta_start, left_coords);

    // solution in polar coordinates
    Vector<double>      u_polar(3, 0.0);
    DenseMatrix<double> grad_u_polar(2,2);
    
    // get the coordinates of this location relative to a coordinate system
    // centred on the left edge of the plate
    Vector<double> relative_coords(2);
    relative_coords[0] = x[0] - left_coords[0];
    relative_coords[1] = x[1] - left_coords[1];

    // QUEHACERES also account for rotation here

    // radius and angle as measured from the edge of the plate
    double r_relative = sqrt(pow(relative_coords[0],2) + pow(relative_coords[1],2));

    // QUEHACERES this will only work for flat plates!!
    // catch slight overshoot of lower plate nodes
    double tol_y = 1.0e-12;
    if ((relative_coords[1] >= 0.0) && (relative_coords[1] < tol_y)
	&& boundary_id == int(Global_Physical_Variables::Inner_plate_boundary))
    {
      // set to small negative number so the angle is ~2pi
      relative_coords[1] = -tol_y;
    }

    if ((relative_coords[1] < 0.0) && (relative_coords[1] > -tol_y)
	&& boundary_id == int(Global_Physical_Variables::Inner_plate_boundary_upper))
    {
      // set to 0 so that the angle is 0
      relative_coords[1] = 0.0;
    }
    
    double theta_relative =
      Additional_Maths_Functions::atan2pi(relative_coords[1], relative_coords[0]);

        
    // get the moffat solution for this relative position in polar coordinates
    moffat_solution(r_relative, theta_relative, A, B, u_polar, grad_u_polar);

    // check if we're at the origin
    bool at_origin = r_relative == 0;
    
    // make sure we've got the right amount of space
    u.resize(3);
    du_dx.resize(2);

    // convert velocity components to Cartesians
    u[0] = u_polar[0] * cos(theta_relative) - u_polar[1] * sin(theta_relative);
    u[1] = u_polar[0] * sin(theta_relative) + u_polar[1] * cos(theta_relative);

    // pressure is a scalar so unchanged by coordinate transform
    u[2] = u_polar[2];
    
    // convert velocity gradient from polar to cartesian
    du_dx = polar_to_cartesian_derivatives_2d(grad_u_polar, u_polar, r_relative, theta_relative);

    // catch the point exactly at the origin
    if(at_origin)
    {
      // from BCs
      u[0] = 0.0;
      u[1] = 0.0; 
      u[2] = -infinity;

      du_dx(0,0) =  infinity;
      du_dx(0,1) = -infinity;
      du_dx(1,0) =  infinity;
      du_dx(1,1) = -infinity;
    }

    // QUEHACERES account for plate velocity
    // u[0] += Global_Physical_Variables::plate_velocity[0];
    // u[1] += Global_Physical_Variables::plate_velocity[1];
    
    // // QUEHACERES debug
    // u[0] = x[0]; // r_relative;
    // u[1] = 0.0;
    // u[2] = 0.0;

    // du_dx(0,0) = 0.0;
    // du_dx(0,1) = 0.0;
    // du_dx(1,0) = 0.0;
    // du_dx(1,1) = 0.0;
    
  }  

  /// \short "Singular" function for right hand edge
  Vector<double> singular_fct_and_gradient_right(const Vector<double>& x,
						 Vector<double>& u, DenseMatrix<double>& du_dx,
						 int boundary_id=-1)
  {
    Vector<double>x_flipped(2);
    
    // flip x coordinate so appears as though the plate is at x>0
    // sort this out at some point, remove left shifting from singular_fct_and_gradient()
    x_flipped[0] = - (x[0]); //  - Global_Physical_Variables::plate_radius);

    // y coordinate unchanged
    x_flipped[1] = x[1];
      
    singular_fct_and_gradient(x_flipped, u, du_dx, boundary_id);

    // account for the flip in the y-axis;
    // u_x, du_x/dy and du_y/dx need to be flipped, but not du_x/dx as the signs cancel.
    u[0] = -u[0];
    
    du_dx(0,1) = -du_dx(0,1);
    du_dx(1,0) = -du_dx(1,0);
    
    return u;
  }
  
  /// \short "Singular" function
  Vector<double> singular_fct(const Vector<double>& x, int boundary_id=-1)
  {
    Vector<double> u(3);
    DenseMatrix<double> du_dx(2,2);
    
    singular_fct_and_gradient(x, u, du_dx, boundary_id);
    
    return u;
  }

  /// \short "Singular" function
  Vector<double> singular_fct_right(const Vector<double>& x, int boundary_id=-1)
  {
    Vector<double> u(3);
    DenseMatrix<double> du_dx(2,2);
    
    singular_fct_and_gradient_right(x, u, du_dx, boundary_id);
    
    return u;
  }
  
  /// \short Gradient of "Singular" function
  DenseMatrix<double> gradient_of_singular_fct(const Vector<double>& x, int boundary_id=-1)
  {
    Vector<double> u(3, 0.0);
    DenseMatrix<double> du_dx(2,2);    
    
    singular_fct_and_gradient(x, u, du_dx, boundary_id);

    return du_dx;
  }

  /// \short Gradient of "Singular" function
  DenseMatrix<double> gradient_of_singular_fct_right(const Vector<double>& x, int boundary_id=-1)
  {
    Vector<double> u(3, 0.0);
    DenseMatrix<double> du_dx(2,2);    
    
    singular_fct_and_gradient_right(x, u, du_dx, boundary_id);

    return du_dx;
  }
  
  /// \short Newtonian stress tensor
  DenseMatrix<double> get_stress(const DenseMatrix<double>& strain_rate,
				 const double& p)
  {
    // \tau_{ij}
    DenseMatrix<double> stress(2,2);

    for(unsigned i=0; i<2; i++)
    {
      for(unsigned j=0; j<2; j++)
      {
	// Newtonian constitutive relation
	stress(i,j) = -p*Additional_Maths_Functions::delta(i,j) + 2.0*strain_rate(i,j);
      }
    }

    return stress;
  }
  
  double get_exact_pressure(const Vector<double>& x)
  {
    // // radius
    // double r = sqrt(x[0]*x[0] + x[1]*x[1]);

    // // polar angle
    // double theta = Additional_Maths_Functions::atan2pi(x[1], x[0]);

    // // analytic result for p
    // double p = -(2.0/sqrt(r)) * ( A * sin(theta/2.0) + 3.0*B*cos(theta/2.0) );

    // // catch infinities and set to a large number
    // if (r == 0)
    // {
    //   // figure out if we're going towards +/- infinity
    //   double sign = Additional_Maths_Functions::sgn( -(A * sin(theta/2.0) + 3.0*B*cos(theta/2.0)) );

    //   // set to a large finite number
    //   p = sign * 1000;
    // }
    
    // return p;
  }

  // computes exact solution at vector position x and returns vector (u,v,p)
  void get_exact_soln_as_vector(const Vector<double>& x, Vector<double>& u)
  {
    // // Radius
    // double r = sqrt(x[0]*x[0] + x[1]*x[1]);

    // double y = x[1];
    
    // // Polar angle
    // double theta = Additional_Maths_Functions::atan2pi(y, x[0]);
    
    // // velocity components for general Panton solution (from Mathematica)
    // double u_r = sqrt(r) * (
    //   -(B*((-3*cos(theta/2.))/2. + (3*cos((3*theta)/2.))/2.))
    //   + A*(sin(theta/2.)/2. - (3*sin((3*theta)/2.))/2.) );
    
    // double u_theta = (-3*sqrt(r)*(A*(-cos(theta/2.) +
    // 				     cos((3*theta)/2.)) - B*(-3*sin(theta/2.)
    // 							     + sin((3*theta)/2.))))/2.;

    // // resize in case an un-initialised vector gets passed in
    // u.resize(3);
    
    // // convert to Cartesians
    // u[0] = u_r * cos(theta) - u_theta * sin(theta);
    // u[1] = u_r * sin(theta) + u_theta * cos(theta);

    // // get the pressure
    // u[2] = get_exact_pressure(x);
  }
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


//==start_of_problem_class============================================
/// Problem class for fibre moving in a 2D box
//====================================================================
template<class ELEMENT> 
class FibreInTwoDBoxProblem : public Problem
{

public:

  /// Constructor
  FibreInTwoDBoxProblem();

  /// Destructor 
  ~FibreInTwoDBoxProblem()
    {
      // hierher: at some point delete things properly and do memory leak
      // check
      // QUEHACERES add back in when we're doing refinable meshes
      // delete Bulk_mesh_pt->spatial_error_estimator_pt();
      delete Bulk_mesh_pt;
    }
 
  /// Update the after solve (empty)
  void actions_after_newton_solve(){}
 
  /// \short Update the problem specs before solve (empty)
  void actions_before_newton_step()
    {
      // QUEHACERES
      
      // residual vector and Jacobian matrix
      DoubleVector r;
      CRDoubleMatrix jac;

      // get 'em
      get_jacobian(r,jac);
  
      ofstream some_file;      
      char filename[100];
      
      sprintf(filename,"%s/jacobian_sparse%i_newton_step_%i.dat", Doc_info.directory().c_str(),
	      Doc_info.number(), Newton_step_counter);      
      some_file.open(filename);
      
      jac.sparse_indexed_output(some_file);

      some_file.close();

      Newton_step_counter++;
    }
  
  // Perform actions after mesh adaptation
  void actions_after_adapt()
    {      
      // Recreate plate nodes      
      duplicate_plate_nodes_and_add_boundary();
      
      // Complete problem setup
      complete_problem_setup();

      // Rebuild global mesh
      rebuild_global_mesh();
  
    }
 
  /// Perform actions after mesh adaptation (empty)
  void actions_before_adapt()
    {
      // Rebuild global mesh
      rebuild_global_mesh();
    }
 
  /// Access function for the specific mesh
  // TriangleMesh<ELEMENT>* mesh_pt()
  RefineableTriangleMesh<ELEMENT>* mesh_pt()  
    {
      return dynamic_cast<RefineableTriangleMesh<ELEMENT>*>(Problem::mesh_pt());
      // return dynamic_cast<TriangleMesh<ELEMENT>*>(Problem::mesh_pt());
    }
 
  /// Doc the solution
  void doc_solution();

  Plate* plate_pt()
    {
      return Plate_pt;
    }

  // function to directly impose the singular amplitude and bypass the 
  // proper calculation
  void impose_fake_singular_amplitude();

  // for debug / validation of singular solution
  void set_values_to_singular_solution();
  
private:

  unsigned Newton_step_counter;
  
  /// hierher Delete face elements and flush meshes
  void delete_face_elements()
  {
    if (CommandLineArgs::command_line_flag_has_been_set
	("--dont_subtract_singularity"))
    {
      return;
    }

    // Loop over the bc elements
    unsigned n_element = Face_mesh_for_bc_pt->nelement();
    for(unsigned e=0; e<n_element; e++)
    {
      // Kill
      delete Face_mesh_for_bc_pt->element_pt(e);
    }
   
    // Wipe the mesh
    Face_mesh_for_bc_pt->flush_element_and_node_storage();

    for(unsigned ising=0; ising<Nsingular_fct; ising++)
    {
      // Loop over the integral face elements
      n_element = Face_mesh_for_singularity_integral_pt[ising]->nelement();
      for(unsigned e=0; e<n_element; e++)
      {
	delete Face_mesh_for_singularity_integral_pt[ising]->element_pt(e);
      }
      Face_mesh_for_singularity_integral_pt[ising]->flush_element_and_node_storage();
    }
  }
    
  /// Create face elements
  void create_face_elements()
  {
    // add Lagrange multipliers to enforce Dirichlet boundary conditions
    if (CommandLineArgs::command_line_flag_has_been_set
	("--enforce_dirichlet_bcs_by_lagrange_multipliers"))
    {            
      // loop over all the outer boundaries
      for(unsigned i_bound=Global_Physical_Variables::First_outer_boundary;
	  i_bound < Global_Physical_Variables::Last_outer_boundary; i_bound++)
      {
	// get the number of elements on this boundary
	unsigned n_element = Bulk_mesh_pt->nboundary_element(i_bound);
          
	// Loop over the bulk elements adjacent to boundary b
	for(unsigned e=0; e<n_element; e++)
	{
	  // Get pointer to the bulk element that is adjacent to boundary b
	  ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
	    Bulk_mesh_pt->boundary_element_pt(i_bound, e));
	      
	  //Find the index of the face of element e along boundary b 
	  int face_index = Bulk_mesh_pt->face_index_at_boundary(i_bound, e);
	  	  
	  // Build the corresponding bc element
	  // N.B. we're calling this without the optional ID, since we have a continuous
	  // Lagrange multiplier field all the way around the outer boundary due to the
	  // continuous no-slip boundary conditions
	  NavierStokesWithSingularityBCFaceElement<ELEMENT>* bc_element_pt =
	    new NavierStokesWithSingularityBCFaceElement<ELEMENT>
	    (bulk_elem_pt, face_index);
            
	  // Tell the element about the singular fct
	  if (!CommandLineArgs::command_line_flag_has_been_set
	      ("--dont_subtract_singularity"))
	  {
	    Vector<ScalableSingularityForNavierStokesElement<ELEMENT>*> sing_el_pt(Nsingular_fct);

	    for(unsigned ising=0; ising<Nsingular_fct; ising++)
	    {
	      sing_el_pt[ising] = dynamic_cast<ScalableSingularityForNavierStokesElement<ELEMENT>*>(
		Singular_fct_element_mesh_pt->element_pt(ising));
	    }
	    
	    bc_element_pt->set_navier_stokes_sing_el_pt(sing_el_pt);
	  }
	  
	  //Add the bc element to the surface mesh
	  Face_mesh_for_bc_pt->add_element_pt(bc_element_pt);
	}
      } // end loop over boundaries
    }

    if (CommandLineArgs::command_line_flag_has_been_set("--dont_subtract_singularity"))
    {
      return;
    }

    // Create the face elements needed to compute the amplitude of
    // the singular function,
   
    // All outer boundaries
    for(unsigned i_bound=Global_Physical_Variables::First_outer_boundary;
	  i_bound < Global_Physical_Variables::Last_outer_boundary; i_bound++)
    {
      unsigned n_element = Bulk_mesh_pt->nboundary_element(i_bound);
      for(unsigned e=0; e<n_element; e++)
      {
	//Create Pointer to bulk element adjacent to the boundary
	ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>
	  (Bulk_mesh_pt->boundary_element_pt(i_bound, e));
       
	//Get Face index of boundary in the bulk element
	int face_index = Bulk_mesh_pt->face_index_at_boundary(i_bound, e);

	for(unsigned ising=0; ising < Nsingular_fct; ising++)
	{
	  //Create corresponding face element
	  NavierStokesWithSingularityBoundaryIntegralFaceElement<ELEMENT>* 
	    boundary_integral_face_element_pt =
	    new NavierStokesWithSingularityBoundaryIntegralFaceElement<ELEMENT>(
	      bulk_elem_pt, face_index, i_bound);

	  // We pass the pointer of singular function element to the face element
	  boundary_integral_face_element_pt->navier_stokes_sing_el_pt() =
	    dynamic_cast<ScalableSingularityForNavierStokesElement<ELEMENT>*>(
	      Singular_fct_element_mesh_pt->element_pt(ising));

	  //Attach it to the mesh
	  Face_mesh_for_singularity_integral_pt[ising]->add_element_pt(boundary_integral_face_element_pt);
	}
      }
    }

    for(unsigned ising=0; ising < Nsingular_fct; ising++)
    {
      // Update the pointer to the face elements (currently needed so
      // this GeneralisedElement can assemble the contributions to the
      // r_C residual from the face elements!
      dynamic_cast<ScalableSingularityForNavierStokesElement<ELEMENT>*>(
	Singular_fct_element_mesh_pt->element_pt(ising))->
	set_mesh_of_face_elements(Face_mesh_for_singularity_integral_pt[ising]);
    }
    
  } // end of create_elements()

  // object containing relevant info for documenting the solution
  DocInfo Doc_info;
  
  // dimensionality of the problem
  unsigned Dim;

  // number of singular functions we're subtracting
  unsigned Nsingular_fct;
  
  /// Do what it says
  void complete_problem_setup();
 
  /// Helper function to apply boundary conditions
  void apply_boundary_conditions();
  
  /// function to duplicate nodes on the plate to allow for jump in pressure
  /// either side of the plate
  void duplicate_plate_nodes_and_add_boundary();

  /// function to attach face elements to the internal boundaries which separate the upper
  /// and lower regions
  void attach_internal_boundary_face_elements();
  
  /// Pointer to the bulk mesh
  RefineableTriangleMesh<ELEMENT>* Bulk_mesh_pt;

  /// Pointer to the "surface" mesh which will have the traction elements
  Mesh* Surface_mesh_pt;

  /// Pointer to the internal surface mesh which will hold face elements along
  /// the internal region boundaries
  Mesh* Internal_boundary_surface_mesh_pt;

  /// Face element mesh for BC (Lagrange multiplier!) 
  Mesh* Face_mesh_for_bc_pt;
  
  /// \short Meshes of face elements used to compute the amplitudes of the singular
  /// functions (one mesh per singular function)
  Vector<Mesh*> Face_mesh_for_singularity_integral_pt;
 
  /// Mesh for (single) element containing singular fct
  Mesh* Singular_fct_element_mesh_pt;
  
  /// Geom object made of mesh
  MeshAsGeomObject* Mesh_as_geom_object_pt;

  // pointer to the plate Geometric Object
  Plate* Plate_pt;
  
}; // end_of_problem_class


//==start_of_constructor==================================================
/// Constructor for FibreInTwoDBoxProblem problem
//========================================================================
template<class ELEMENT>
FibreInTwoDBoxProblem<ELEMENT> :: FibreInTwoDBoxProblem() : Newton_step_counter(0)
{
  // set the dimensionality
  Dim = 2;

  // number of singular functions that we're subtracting
  Nsingular_fct = 2;
  
  // Set up doc info  
  // ---------------

  // don't run if the output folder doesn't exist
  Doc_info.enable_error_if_directory_does_not_exist();
  
  // Set output directory
  Doc_info.set_directory(Global_Physical_Variables::dir.c_str());
 
  // Step number
  Doc_info.number() = 0;

  if (CommandLineArgs::command_line_flag_has_been_set("--circular_plate"))
  {
    // shorthand
    double r_c = Global_Physical_Variables::plate_radius_of_curvature;
  
    // centre of the circle; central along x-axis of plate, and with
    // a y coordinate such that the two plate ends touch the edge of the circle
    double x_c = 0.0;
    double y_c = sqrt( pow(r_c,2) - pow(Global_Physical_Variables::plate_radius, 2) );

    // intrinsic coordinates of the start and end of the plate
    double s_plate_start = atan2(-y_c, -Global_Physical_Variables::plate_radius);
    double s_plate_end   = atan2(-y_c, +Global_Physical_Variables::plate_radius);
  
    // Make GeomObject representing a circular plate
      Plate* plate_pt =
      new CircularPlate(x_c, y_c, r_c, s_plate_start, s_plate_end,
			Global_Physical_Variables::high_res_region_zeta);
  }
  else 
  {
    // default is flat plate
    
    Vector<double> start_coords(2, 0.0);
    Vector<double> end_coords(2, 0.0);
  
    start_coords[0] = - Global_Physical_Variables::plate_radius;
    end_coords[0]   = + Global_Physical_Variables::plate_radius;

  
    Plate_pt = new FlatPlate(start_coords, end_coords,
			     Global_Physical_Variables::high_res_region_zeta);
  }

  
  // Build the mesh
  Bulk_mesh_pt = Global_Physical_Variables::build_the_mesh<ELEMENT>
    (Global_Physical_Variables::element_area, Plate_pt);

  // Let's have a look at the boundary enumeration
  Bulk_mesh_pt->output_boundaries("boundaries.dat");

  // Set error estimator for bulk mesh
  Z2ErrorEstimator* error_estimator_pt = new Z2ErrorEstimator;
  Bulk_mesh_pt->spatial_error_estimator_pt() = error_estimator_pt;

  // Set element size limits
  Bulk_mesh_pt->max_element_size()=0.1;
  Bulk_mesh_pt->min_element_size()=1e-30;
  Bulk_mesh_pt->max_permitted_error()=0.005;
  Bulk_mesh_pt->min_permitted_error()=0.0;

  // duplicate nodes above the plate to allow for a pressure jump across the plate
  duplicate_plate_nodes_and_add_boundary();

  // Add sub-meshes
  add_sub_mesh(Bulk_mesh_pt);

  // create the internal surface mesh which will allow the solution to be computed along the
  // plate and along the boundaries which divide the upper and lower regions
  Internal_boundary_surface_mesh_pt = new Mesh;
  
  // attach face elements to bulk elements in the upper region on the internal boundaries
  attach_internal_boundary_face_elements();
  
  add_sub_mesh(Internal_boundary_surface_mesh_pt);

  // check we're not doing pure FE
  if (!CommandLineArgs::command_line_flag_has_been_set("--dont_subtract_singularity"))
  {

    // ================================================================
    // QUEHACERES need to add additional sinuglarities here...
    // for the time being we'll do just one edge for the case 
    
    // Create element that stores the singular fct and its amplitude
    //---------------------------------------------------------------
    ScalableSingularityForNavierStokesElement<ELEMENT>* el_pt =
      new ScalableSingularityForNavierStokesElement<ELEMENT>;

    ScalableSingularityForNavierStokesElement<ELEMENT>* right_el_pt =
      new ScalableSingularityForNavierStokesElement<ELEMENT>;
    
    // Pass fct pointers:
    el_pt->unscaled_singular_fct_pt() = &Analytic_Functions::singular_fct;
    el_pt->gradient_of_unscaled_singular_fct_pt() =
      &Analytic_Functions::gradient_of_singular_fct;

    right_el_pt->unscaled_singular_fct_pt() = &Analytic_Functions::singular_fct_right;
    right_el_pt->gradient_of_unscaled_singular_fct_pt() =
      &Analytic_Functions::gradient_of_singular_fct_right;
    
    // Add to mesh
    Singular_fct_element_mesh_pt = new Mesh;
    Singular_fct_element_mesh_pt->add_element_pt(el_pt);

    // QUEHACERES don't add the right one yet
    Singular_fct_element_mesh_pt->add_element_pt(right_el_pt);
    
    add_sub_mesh(Singular_fct_element_mesh_pt);

    // Create face elements that compute contribution to amplitude residual
    //---------------------------------------------------------------------
    for(unsigned ising=0; ising < Nsingular_fct; ising++)
    {
      Face_mesh_for_singularity_integral_pt.push_back(new Mesh);
    }
  }

  // Create face elements for imposition of BC
  Face_mesh_for_bc_pt = new Mesh;

  // Build the face elements
  create_face_elements();

  // Add 'em to mesh
  add_sub_mesh(Face_mesh_for_bc_pt);
  
  // Build global mesh
  build_global_mesh();
  
  // Complete problem setup (apply boundary conditions and set pointers)
  complete_problem_setup();

  oomph_info << "number of elements on plate: "
	     << Bulk_mesh_pt->nboundary_element(Global_Physical_Variables::Inner_plate_boundary)
	     << "\nnumber of nodes on plate:    "
	     << Bulk_mesh_pt->nboundary_node(Global_Physical_Variables::Inner_plate_boundary)
	     << "\n\n";
  
  // Setup equation numbering scheme
  oomph_info <<"Number of equations: " 
             << this->assign_eqn_numbers() 
             << std::endl;

} // end_of_constructor


template <class ELEMENT>
class InternalBoundarySurfaceElement :  public virtual FaceGeometry<ELEMENT>,
					public virtual FaceElement
{
public:

  InternalBoundarySurfaceElement(FiniteElement* const &element_pt,
				 const int &face_index) :    
    FaceGeometry<ELEMENT>()
    {
      element_pt -> build_face_element(face_index, this);
    }

  /// \short Output function
  void output(std::ostream &outfile, const unsigned &n_plot)
  {    
    // Get pointer to assocated bulk element
    ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(bulk_element_pt());
   
    // Elemental dimension
    unsigned dim_el=dim();
   
    //Local coordinates
    Vector<double> s(dim_el);

    Vector<double> s_bulk(dim_el+1);
       
    //Calculate the Eulerian coordinates
    Vector<double> x(dim_el+1, 0.0);

    // Velocity from bulk element
    Vector<double> veloc(dim_el+1);

    // pressure
    double p;
    
    // // Tecplot header info
    // outfile << this->tecplot_zone_string(n_plot);
   
    // get number of plot points
    unsigned num_plot_points = this->nplot_points(n_plot);

    // Loop over plot points
    for (unsigned iplot=0; iplot<num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      this->get_s_plot(iplot,n_plot,s);
     
      this->get_local_coordinate_in_bulk(s, s_bulk);
     
      //Get x position from bulk
      bulk_el_pt->interpolated_x(s_bulk, x);

      // get velocity
      bulk_el_pt->interpolated_u_nst(s_bulk, veloc);

      // get pressure
      p = bulk_el_pt->interpolated_p_nst(s_bulk);

      // output the coordinates
      for(unsigned i=0; i<dim_el+1; i++)
      {
	outfile << x[i] << "  ";
      }

      // output the velocity
      for(unsigned i=0; i<dim_el+1; i++)
      {
	outfile << veloc[i] << " ";
      }

      // output pressure
      outfile << p;

      // done with this plot point
      outfile << std::endl;
    }

    // this->write_tecplot_zone_footer(outfile,n_plot);    
  }

private:
  
  /// Add the element's contribution to its residual vector
   inline void fill_in_contribution_to_residuals(Vector<double> &residuals)
   {
   }

   /// \short Add the element's contribution to its residual vector and its
   /// Jacobian matrix
   inline void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                                DenseMatrix<double> &jacobian)
   {
   }
};

//============start_of_attach_internal_boundary_face_elements============
/// Attach face elements to the bulk elements which are in the lower region
/// and have one side on one of the internal boundaries (either region boundaries
/// or the upper plate boundary)
//=======================================================================
template<class ELEMENT>
void FibreInTwoDBoxProblem<ELEMENT> :: attach_internal_boundary_face_elements()
{
  // vector of regions in lower half of domain
  Vector<unsigned> lower_regions(3);
  lower_regions[0] = Global_Physical_Variables::Region_hi_res_left_lower;
  lower_regions[1] = 0; // the lower region
  lower_regions[2] = Global_Physical_Variables::Region_hi_res_right_lower; 
      
  // iterator to the vector of lower regions
  Vector<unsigned>::iterator region_pt; 

  unsigned end_boundary = Global_Physical_Variables::Inner_plate_boundary;

  end_boundary = Global_Physical_Variables::Inner_region_boundary_right;
  
  // loop over the boundaries dividing the top/bottom of the domain
  for(unsigned b = Global_Physical_Variables::Inner_region_boundary_left;
      b <= end_boundary; b++ )
  {
    for (region_pt = lower_regions.begin(); region_pt < lower_regions.end(); region_pt++)
    {    
      // How many bulk fluid elements are adjacent to this boundary b in the upper region?
      unsigned n_element = Bulk_mesh_pt->nboundary_element_in_region(b, *region_pt);

      // loop over them
      for(unsigned e = 0; e < n_element; e++)
      {
	// Get pointer to the bulk fluid element that is 
	// adjacent to boundary b in region 0
	ELEMENT* bulk_elem_pt =
	  dynamic_cast<ELEMENT*>(Bulk_mesh_pt->boundary_element_in_region_pt(b, *region_pt, e));
       
	//Find the index of the face of element e along boundary b in region 0
	int face_index = Bulk_mesh_pt->
	  face_index_at_boundary_in_region(b, *region_pt, e);

	// make a new face element
	InternalBoundarySurfaceElement<ELEMENT>* surface_element_pt =
	  new InternalBoundarySurfaceElement<ELEMENT>(bulk_elem_pt, face_index);

	// Add it to the surface sub mesh
	Internal_boundary_surface_mesh_pt->add_element_pt(surface_element_pt);
      }
    }
  }
}

template<class ELEMENT>
void FibreInTwoDBoxProblem<ELEMENT> :: duplicate_plate_nodes_and_add_boundary()
{
  
  // Want elements in the upper region on the plate
  unsigned b         = Global_Physical_Variables::Inner_plate_boundary;
  unsigned region_id = Global_Physical_Variables::Region_upper;
  unsigned nel_bulk  = Bulk_mesh_pt->nboundary_element_in_region(b,region_id);

  // we also need to duplicate the nodes which are in the high-res regions above the plate,
  unsigned nel_high_res_left =
    Bulk_mesh_pt->nboundary_element_in_region( b,
					       Global_Physical_Variables::Region_hi_res_left_upper);
  unsigned nel_high_res_right =
    Bulk_mesh_pt->nboundary_element_in_region( b,
					       Global_Physical_Variables::Region_hi_res_right_upper);


  // total number of elements above the plate
  unsigned nel_total = nel_bulk + nel_high_res_left + nel_high_res_right;
  
   // ====================================================
#ifdef DEBUG
  {
    unsigned nnode_total = Bulk_mesh_pt->nnode();
      	

    unsigned n_plate_nodes = 0;

    oomph_info << "\n-----------------------------\n"
  	       << "Before duplication, nnode = " << nnode_total << "\n"
  	       << "\n-----------------------------\n\n";

    for(unsigned i=0; i< nnode_total; i++)
    {
      Node* node_pt = Bulk_mesh_pt->node_pt(i);
      oomph_info << "node "<< i << ": (" << node_pt->x(0) << ", " << node_pt->x(1) << ")";

      if( node_pt->is_on_boundary(b) )
      {	      
  	n_plate_nodes++;

  	oomph_info << " [Plate] \n";
      }
      else
      {
  	oomph_info << "\n";
      }
    }

    oomph_info << "\n number of nodes on plate: " << n_plate_nodes << "\n\n";
  }
#endif
  
  // ====================================================
  
  oomph_info << "Number of boundaries before resize: " << Bulk_mesh_pt->nboundary() << std::endl;
      
  // increase the number of boundaries in the mesh to accomodate the extra plate boundary
  Bulk_mesh_pt->set_nboundary(Bulk_mesh_pt->nboundary()+1);

  oomph_info << "Number of boundaries after resize:  " << Bulk_mesh_pt->nboundary() << "\n\n";

  oomph_info << "Total number of nodes in mesh before duplication: "
	     << Bulk_mesh_pt->nnode() << "\n\n" << std::flush; 
  
  // map to keep track of nodes we've already duplicated;
  // map is original node -> new node
  std::map<Node*,Node*> existing_duplicate_node_pt;

  // loop over the elements above the plate
  for (unsigned e=0; e<nel_total; e++)
  {
    // we need to fiddle the element number here, as we're traversing 3 different regions
    // (left high-res, the upper bulk, and the right high-res regions) so we need a numbering
    // which accounts for the region we're in, and need to keep track of the regions
    unsigned e_region;   
    
    if( e < nel_high_res_left )
    {
      // we're moving to the right, so if we're at the start in the left high-res region,
      // the loop variable is correct, nothing else to do
      e_region = e;
      region_id = Global_Physical_Variables::Region_hi_res_left_upper;
    }
    else if( e >= nel_high_res_left &&
	     e <  nel_high_res_left + nel_bulk )
    {
      // we're in the bulk region now, so we subtract off the number
      // of elements in the left region, so that we start again
      // counting from the zeroth element in the bulk region
      e_region = e - nel_high_res_left;
      region_id = Global_Physical_Variables::Region_upper;
    }
    else
    {
      // we're in the right high-res region now, so we subtract off the
      // number of elements in both the left high-res region and bulk region
      e_region = e - (nel_high_res_left + nel_bulk);
      region_id = Global_Physical_Variables::Region_hi_res_right_upper;
    }
    
    // get a pointer to this element
    FiniteElement* el_pt = Bulk_mesh_pt->boundary_element_in_region_pt(b, region_id, e_region);

    // get the number of nodes in this element
    unsigned nnode = el_pt->nnode();

    // loop over nodes
    for(unsigned j=0; j<nnode; j++)
    {
      // get a pointer to the node
      Node* nod_pt = el_pt->node_pt(j);

      // QUEHACERES probably delete this, we known it's on the boundary because of the function we
      // used to get the element pointer above
      // check if this node is on the plate
      bool is_plate_node = nod_pt->is_on_boundary(b);

      // QUEHACERES debug
      Vector<double> x(2);
      x[0] = nod_pt->x(0);
      x[1] = nod_pt->x(1);
	
      // QUEHACERES we're back to not duplicating these for the time being
      // // QUEHACERES for now we will duplicate the edge nodes too to allow for +/- infinity as
      // // the singularity is approached from above and below
      if( is_plate_node &&
	  !nod_pt->is_on_boundary(Global_Physical_Variables::Inner_region_boundary_left)
	  &&
	  !nod_pt->is_on_boundary(Global_Physical_Variables::Inner_region_boundary_right)
	)	
      {	    
	// Find this original node in the map; if we find it
	// it's already been duplicated earlier; in that case
	// use that duplicate rather than creating another one.
	std::map<Node*,Node*>::iterator it =
	  existing_duplicate_node_pt.find(nod_pt);
	    
	unsigned n_dim           = nod_pt->ndim();
	unsigned n_position_type = nod_pt->nposition_type();
	unsigned n_value         = nod_pt->nvalue();

	// check if we've already duplicated this node
	if (it == existing_duplicate_node_pt.end())
	{	    
	  // create a new node
	  el_pt->node_pt(j) = new BoundaryNode<Node>(n_dim, n_position_type, n_value);
         
	  // It has the same coordinate; hierher add history values too
	  // when implementing as Navier Stokes
	  for (unsigned i=0; i<n_dim; i++)
	  {
	    el_pt->node_pt(j)->x(i) = nod_pt->x(i);
	  }
         
	  // ...and the same values
	  for (unsigned i=0; i<n_value; i++)
	  {
	    el_pt->node_pt(j)->set_value(i, nod_pt->value(i));
	  }

	  // =================

	  // ID of new boundary
	  unsigned new_boundary_id = Global_Physical_Variables::Inner_plate_boundary_upper;
	      
	  // add new node to the new boundary
	  // el_pt->node_pt(j)->add_to_boundary(new_boundary_id);
	  // calling this both calls the node pointers function to tell it it's on the boundary,
	  // and also updates the mesh's list of boundary nodes
	  Bulk_mesh_pt->add_boundary_node( new_boundary_id, el_pt->node_pt(j) );
	      
	  // Get/set boundary coordinates
	  if ( nod_pt->boundary_coordinates_have_been_set_up() )
	  {
	    // get number of coordinates on the original plate boundary
	    unsigned ncoords = nod_pt->ncoordinates_on_boundary(b);
		
	    Vector<double> boundary_zeta(ncoords);

	    // get 'em from original plate boundary
	    nod_pt->get_coordinates_on_boundary(b, boundary_zeta);
		
	    // set 'em for new plate boundary
	    el_pt->node_pt(j)->set_coordinates_on_boundary(new_boundary_id, boundary_zeta);

	    // shorthand
	    unsigned obu = Global_Physical_Variables::Outer_boundary_upper;

	    // also want to add new node to the outer boundary if the original one was there too
	    if(nod_pt->is_on_boundary(obu))
	    {
	      // add new node to the upper outer boundary
	      el_pt->node_pt(j)->add_to_boundary(obu);
	      
	      // get number coordinates on outer boundary
	      ncoords = nod_pt->ncoordinates_on_boundary(obu);

	      // get 'em
	      nod_pt->get_coordinates_on_boundary(obu, boundary_zeta);

	      // set 'em
	      el_pt->node_pt(j)->set_coordinates_on_boundary(obu, boundary_zeta);
	    }
	  }
	  else
	  {
	    // hierher throw? (Doesn't happen at the moment, i.e. 
	    // when this diagnostic was finally commented out)
             
	    oomph_info << "No boundary coordinates have been set up"
		       << " for new local node " << j
		       << " at : "
		       << nod_pt->x(0) << " "
		       << nod_pt->x(1)
		       << std::endl;
	  }	    
	  
	  // add  new node to the bulk mesh
	  Bulk_mesh_pt->add_node_pt(el_pt->node_pt(j));
	      
	  // Keep track, add entry mapping from the old node to the new node
	  existing_duplicate_node_pt[nod_pt] = el_pt->node_pt(j);
	}
      }
    }
  }
  
  // QUEHACERES probably need this since we've fiddled the nodes on the plate boundary
  Bulk_mesh_pt->setup_boundary_element_info();

  // ================================================================================================
  // Now we've fixed the boundary elements, i.e. the ones who have a edge along the boundary,
  // we need to fix the adjacent elements who only have a single vertex on the boundary (which
  // are not considered boundary elements so don't get returned by boundary_element_in_region_pt(...) )
	
  // get the number of elements in the upper region	
  unsigned nel_upper_bulk_region =
    Bulk_mesh_pt->nregion_element(Global_Physical_Variables::Region_upper);

  // number of elements in the upper part of the left high-res region
  unsigned nel_upper_left_high_res_region =
    Bulk_mesh_pt->nregion_element(Global_Physical_Variables::Region_hi_res_left_upper);

  // number of elements in the upper part of the right high-res region
  unsigned nel_upper_right_high_res_region =
    Bulk_mesh_pt->nregion_element(Global_Physical_Variables::Region_hi_res_right_upper);

  // total number of elements in the upper half of the domain
  unsigned nel_total_upper_half = nel_upper_bulk_region +
    nel_upper_left_high_res_region + nel_upper_right_high_res_region;
    
  for(unsigned e=0; e<nel_total_upper_half; e++)
  {
    // same trick as above, although we're now grabbing elements in a region, not
    // just the ones that are in a region *and* on the plate boundary, so we need to
    // adjust the region ID here as well as the element counter
    unsigned e_region;

    if( e < nel_upper_left_high_res_region )
    {
      // we're moving to the right, so if we're at the start in the left high-res region,
      // the loop variable is correct, nothing else to do
      e_region = e;
     region_id = Global_Physical_Variables::Region_hi_res_left_upper;
    }
    else if( e >= nel_upper_left_high_res_region &&
	     e <  nel_upper_left_high_res_region + nel_upper_bulk_region )
    {
      // we're in the bulk region now, so we subtract off the number
      // of elements in the left region, so that we start again
      // counting from the zeroth element in the bulk region
      e_region = e - nel_upper_left_high_res_region;
      region_id = Global_Physical_Variables::Region_upper;
    }
    else
    {
      // we're in the right high-res region now, so we subtract off the
      // number of elements in both the left high-res region and bulk region
      e_region = e - (nel_upper_left_high_res_region + nel_upper_bulk_region);
      region_id = Global_Physical_Variables::Region_hi_res_right_upper;      
    }
        
    // get a pointer to the element
    FiniteElement* el_pt = Bulk_mesh_pt->region_element_pt(region_id, e_region);

    // get number of nodes
    unsigned nnode = el_pt->nnode();

    // now loop over its nodes and check if any are on the boundary
    for(unsigned j=0; j<nnode; j++)
    {
      // get a pointer to this node
      Node* node_pt = el_pt->node_pt(j);

      // get coordinates of this node
      unsigned x_coord = node_pt->x(0);
      unsigned y_coord = node_pt->x(1);
      
      
      // check if it's on the original boundary (the one's we're interested in won't be on the
      // new plate boundary yet)
      if( node_pt->is_on_boundary(Global_Physical_Variables::Inner_plate_boundary) &&
         !node_pt->is_on_boundary(Global_Physical_Variables::Inner_region_boundary_left)
	  &&
	 !node_pt->is_on_boundary(Global_Physical_Variables::Inner_region_boundary_right)
	)
      {	      
	// check our map to find out the correct new node	      
	std::map<Node*,Node*>::iterator it =
	  existing_duplicate_node_pt.find(node_pt);

	// did we find it?
	if (it != existing_duplicate_node_pt.end())
	{
	  // update the node pointer
	  el_pt->node_pt(j) = existing_duplicate_node_pt[node_pt];

	  // should be all we need to do, the boundary info, coordinates, values etc. have already
	  // been sorted out when we created the new node
	}	
	else
	{
	  // if we don't find it something weird has happened, because we have a node
	  // on the original plate boundary that didn't get duplicated
	  oomph_info << "\nNode: (" << x_coord << ", " << y_coord << "); "
		     << "Weird tingz have happened here; I've found a node which \n"
		     <<" is apparently on the original plate boundary, but \n"
		     << "which didn't get duplicated\n\n";
	}
      }
    }    
  }
  
  // QUEHACERES debug, final check
  oomph_info << "Total number of nodes in mesh after duplication: "
	     << Bulk_mesh_pt->nnode() << "\n\n"; 
      
} // end duplicate_plate_nodes_and_add_boundary()

//==start_of_complete======================================================
/// Set boundary condition, and complete the build of
/// all elements
//========================================================================
template<class ELEMENT>
void FibreInTwoDBoxProblem<ELEMENT>::complete_problem_setup()
{    
  // Loop over the elements to set up element-specific
  // things that cannot be handled by constructor
  unsigned n_el=Bulk_mesh_pt->nelement();
  for (unsigned e=0; e<n_el; e++)
  {
    ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(
      Bulk_mesh_pt->element_pt(e));

    // Tell bulk element which function computes the stress
    bulk_el_pt->stress_fct_pt() = &Analytic_Functions::newtonian_stress;

    if (!CommandLineArgs::command_line_flag_has_been_set("--dont_subtract_singularity"))
    {
      // get the number of singular functions to subtract
      unsigned nsingular_fct = Singular_fct_element_mesh_pt->nelement();
      for(unsigned ising=0; ising<Nsingular_fct; ising++)
      {
	// Tell the bulk element about each singular fct
	bulk_el_pt->add_singular_fct_pt(
	  dynamic_cast<TemplateFreeScalableSingularityForNavierStokesElement*>(
	    Singular_fct_element_mesh_pt->element_pt(ising)) );
      }
    }
  }
  
  if (!CommandLineArgs::command_line_flag_has_been_set("--dont_subtract_singularity"))
  { 
    //BC elements
    unsigned n_element = Face_mesh_for_bc_pt->nelement();
    for(unsigned e=0; e<n_element; e++)
    {
      // Upcast from GeneralisedElement to the present element
      NavierStokesWithSingularityBCFaceElement<ELEMENT>* el_pt =
	dynamic_cast<NavierStokesWithSingularityBCFaceElement<ELEMENT>*>(
	  Face_mesh_for_bc_pt->element_pt(e));

      Vector<ScalableSingularityForNavierStokesElement<ELEMENT>*> sing_fct(Nsingular_fct);

      for(unsigned ising=0; ising<Nsingular_fct; ising++)
      {
	sing_fct[ising] = dynamic_cast<ScalableSingularityForNavierStokesElement<ELEMENT>*>(
	  Singular_fct_element_mesh_pt->element_pt(ising));
      }
      
      // Tell the element about the singular functions
      el_pt->set_navier_stokes_sing_el_pt(sing_fct);
    }
  }
  
  // Apply bcs
  apply_boundary_conditions();

  // Complete the build of all elements so they are fully functional

  //Find number of elements in mesh
  unsigned n_element = Bulk_mesh_pt->nelement();

  // Loop over the elements to set up element-specific 
  // things that cannot be handled by constructor
  for(unsigned e=0; e<n_element; e++)
  {
    // Upcast from GeneralisedElement to the present element
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

    //Set the Reynolds number
    el_pt->re_pt() = &Global_Physical_Variables::Re;
  } // end loop over elements
}

//==start_of_apply_bc=====================================================
/// Helper function to apply boundary conditions
//========================================================================
template<class ELEMENT>
void FibreInTwoDBoxProblem<ELEMENT>::apply_boundary_conditions()
{
  // local shorthand for linear and angular velocites
  double vel_x          = Global_Physical_Variables::plate_velocity[0];
  double vel_y          = Global_Physical_Variables::plate_velocity[1];
  double boundary_vel_x = Global_Physical_Variables::boundary_velocity[0];
  double boundary_vel_y = Global_Physical_Variables::boundary_velocity[1];
	
  double omega_plate    = Global_Physical_Variables::angular_velocity_plate;
  double omega_boundary = Global_Physical_Variables::angular_velocity_boundary;
  
  // Set the boundary conditions for this problem: All nodes are
  // free by default -- just pin the ones that have Dirichlet conditions
  // here.
  unsigned num_bound = Bulk_mesh_pt->nboundary();
  for(unsigned ibound=0; ibound<num_bound; ibound++)
  {       
      unsigned num_nod = Bulk_mesh_pt->nboundary_node(ibound);
      for (unsigned inod=0; inod<num_nod; inod++)
      {
	// get pointer to this boundary node
	Node* node_pt = Bulk_mesh_pt->boundary_node_pt(ibound,inod);

	// get its location
	Vector<double> x(2);	
	x[0] = node_pt->x(0);
	x[1] = node_pt->x(1);
		
	// compute distance from centre of rotation to this node
	double r_x   = x[0] - Global_Physical_Variables::centre_of_rotation[0];
	double r_y   = x[1] - Global_Physical_Variables::centre_of_rotation[1];
	
	// set BCs on the plate boundaries
	if ( ibound == Global_Physical_Variables::Inner_plate_boundary ||
	     ibound == Global_Physical_Variables::Inner_plate_boundary_upper )
	{
	  // velocity of the node is \bm u + \bm\Omega\times\bm r,
	  // where \bm r is measured from the centre of rotation
	  double u = vel_x - omega_plate * r_y;
	  double v = vel_y + omega_plate * r_x;
	  
	  // pin horizontal and vertical velocities on the plate
	  node_pt->pin(0);
	  node_pt->pin(1);
	
	  // set 'em // QUEHACERES
	  node_pt->set_value(0, u);
	  node_pt->set_value(1, v);
	}
	// set BCs on the outer boundaries
	else if( ibound == Global_Physical_Variables::Outer_boundary_upper ||
		 ibound == Global_Physical_Variables::Outer_boundary_lower )
	{
	  // pin the boundary velocities
	  node_pt->pin(0);
	  node_pt->pin(1);
	  
	  // no slip and no penetration on outer walls
	  node_pt->set_value(0, boundary_vel_x);
	  node_pt->set_value(1, boundary_vel_y);
	}
	else 
	{
	  break;
	}
      }
      
  } // end loop over boundaries
  
  // Free values on boundary if Lagrange multiplier is used
  // to enforce bc for sum of FE soln and singuar fct
  if (CommandLineArgs::command_line_flag_has_been_set("--enforce_dirichlet_bcs_by_lagrange_multipliers"))
  {
    // Now unpin nodal values where the bc conditions are enforced
    // by Lagrange multipliers to ensure that the sum of FE and singular
    // solution is correct
    unsigned nel = Face_mesh_for_bc_pt->nelement();
    for (unsigned e=0; e<nel; e++)
    {
      // Get element
      NavierStokesWithSingularityBCFaceElement<ELEMENT>* el_pt =
	dynamic_cast<NavierStokesWithSingularityBCFaceElement<ELEMENT>*>(
	  Face_mesh_for_bc_pt->element_pt(e));
     
      // Specify desired nodal values for compound solution
      unsigned nnod = el_pt->nnode();

      // matrix to store velocities at each boundary node
      DenseMatrix<double> nodal_boundary_value(nnod, Dim);

      // Unpin the FE part of the solution
      for (unsigned j=0; j<nnod; j++)
      {
	for(unsigned i=0; i<Dim; i++)
	{
	  el_pt->unpin_u_fe_at_specified_local_node(j, i);
	}
		
	Node* node_pt = el_pt->node_pt(j);

	Vector<double> x(2);
	x[0] = node_pt->x(0);
	x[1] = node_pt->x(1);

	// find the boundaries this node is on
	std::set<unsigned>* boundaries_set;
	node_pt->get_boundaries_pt(boundaries_set);

	// assuming BCs are continuous so doesn't matter if node is on multiple
	// boundaries, just take the first
	std::set<unsigned>::iterator boundaries_iter = (*boundaries_set).begin();
	unsigned first_boundary = *boundaries_iter;
	
	// get velocity from boundary conditions at this point
	Vector<double> u(Dim);
	// u = Global_Physical_Variables::u_BC(first_boundary);
	u[0] = boundary_vel_x;
	u[1] = boundary_vel_y;
	
	// assign to the matrix of nodal values
	for(unsigned i=0; i<Dim; i++)
	{
	  nodal_boundary_value(j,i) = u[i];
	}
      }
      
      // Tell the element about these nodal boundary values
      el_pt->set_nodal_boundary_values(nodal_boundary_value);
      
    } // end loop over BC face elements
    
  } // end if(enforced-by-lagrange-mpy)
 
} // end set bc


//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void FibreInTwoDBoxProblem<ELEMENT>::doc_solution()
{
  ofstream some_file;
  char filename[100];

  // Number of plot points
  unsigned npts = 2;

  // save memory by only outputting the hi-res interpolation if requested
  if (CommandLineArgs::command_line_flag_has_been_set("-output_hi_res_soln") )
  {
    npts = 10;
    
    // Output solution
    sprintf(filename,"%s/soln%i.dat", Doc_info.directory().c_str(),
	    Doc_info.number());
    some_file.open(filename);
    Bulk_mesh_pt->output(some_file,npts);
    some_file.close();
  }

  // Output solution just using vertices so we can see the mesh
  sprintf(filename,"%s/coarse_soln%i.dat", Doc_info.directory().c_str(),
	  Doc_info.number());
  some_file.open(filename);
  
  npts = 2;
  Bulk_mesh_pt->output(some_file,npts);
  some_file.close();
  
  // Plot "extended solution" showing contributions; also work out
  // average element size
  double av_el_size = 0.0;
  sprintf(filename,"%s/extended_soln%i.dat",Doc_info.directory().c_str(),
          Doc_info.number());
  some_file.open(filename);
  unsigned nel=Bulk_mesh_pt->nelement();
  for (unsigned e=0; e<nel; e++)
  {
    npts=2;
    ELEMENT* el_pt= dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
    
    el_pt->output_with_various_contributions(some_file, npts);
    av_el_size += el_pt->size();
  }
  some_file.close();
  av_el_size/=double(nel);
  
  // Get error
  double error, norm; 
  sprintf(filename,"%s/error%i.dat", Doc_info.directory().c_str(),
          Doc_info.number());
  some_file.open(filename);
  Bulk_mesh_pt->compute_error(some_file,
                              Analytic_Functions::get_exact_soln_as_vector,
                              error, norm);
  some_file.close();

  // Doc error norm:
  oomph_info << "\n av el size, av h, Ndof, # bulk els, Norm of error    : "   
             << av_el_size << " " 
             << sqrt(av_el_size) << " " 
             << ndof() << " " 
             << Bulk_mesh_pt->nelement() << " " 
             << sqrt(error) << std::endl;
  oomph_info << "Norm of solution : " << sqrt(norm) << std::endl << std::endl;
  oomph_info << std::endl;

  // Exact solution
  sprintf(filename,"%s/exact_soln%i.dat", Doc_info.directory().c_str(),
          Doc_info.number());
  some_file.open(filename);
  
  unsigned nplot = 2;
  Bulk_mesh_pt->output_fct(some_file, nplot,
                           Analytic_Functions::get_exact_soln_as_vector);
  some_file.close();

  // ==============================
  // Plot solution along two axes

  // How many surface elements are in the surface mesh
  unsigned n_element =  Internal_boundary_surface_mesh_pt->nelement();
  
  // Loop over the surface elements
  for(unsigned e = 0; e < n_element; e++)
  {
    // FaceElement* surface_element_pt = 
    //   dynamic_cast<FaceElement*>(Internal_boundary_surface_mesh_pt->element_pt(e));

    // QUEHACERES do something useful here!
  }
  
  // open x-axis file  
  sprintf(filename,"%s/soln_along_x_axis%i.dat", Doc_info.directory().c_str(),
	  Doc_info.number());
  
  some_file.open(filename);

  // output
  nplot = 2;
  Internal_boundary_surface_mesh_pt->output(some_file, nplot);
  some_file.close();

  // residual vector and Jacobian matrix
  DoubleVector r;
  CRDoubleMatrix jac;

  // get 'em
  get_jacobian(r,jac);

  if (CommandLineArgs::command_line_flag_has_been_set("--output_jacobian_sparse"))
  {
    sprintf(filename,"%s/jacobian_sparse%i.dat", Doc_info.directory().c_str(),
	    Doc_info.number());      
    some_file.open(filename);
      
    jac.sparse_indexed_output(some_file);

    some_file.close();
  }

  if (CommandLineArgs::command_line_flag_has_been_set("--describe_dofs"))
  {
    sprintf(filename,"%s/describe_dofs.dat", Doc_info.directory().c_str());
    some_file.open(filename);
    some_file << "test\n";
    describe_dofs(some_file);
    
    some_file.close();
  }
  if (CommandLineArgs::command_line_flag_has_been_set("--describe_nodes"))
  {
    sprintf(filename,"%s/describe_nodes.dat", Doc_info.directory().c_str());
    some_file.open(filename);

    for(unsigned j=0; j<Bulk_mesh_pt->nnode(); j++)
    {
      // grab the node
      Node* node_pt = Bulk_mesh_pt->node_pt(j);

      // get it's coordinates
      double x = node_pt->x(0);
      double y = node_pt->x(1);

      some_file << j << " " << x << " " << y << std::endl;
    }
    
    some_file.close();
  }

  if (!CommandLineArgs::command_line_flag_has_been_set
      ("--dont_subtract_singularity"))
  {
    for(unsigned ising=0; ising<Nsingular_fct; ising++)
    {
      oomph_info 
	<< "Amplitude of singular function " << ising << ": "
	<< dynamic_cast<TemplateFreeScalableSingularityForNavierStokesElement*>(
	  Singular_fct_element_mesh_pt->element_pt(ising))-> 
	amplitude_of_singular_fct() << std::endl;
    }
  }

  if (CommandLineArgs::command_line_flag_has_been_set("--output_z2_error"))
  {
    // grab the error in each element from the Z2 estimator
    Vector<double> elemental_error(Bulk_mesh_pt->nelement());
    Mesh* mesh_pt = dynamic_cast<Mesh*>(Bulk_mesh_pt);
    
    Bulk_mesh_pt->spatial_error_estimator_pt()->get_element_errors(mesh_pt, elemental_error);

    double z2_global_sum = 0;
    
    // sum the errors to get a global measure    
    for(Vector<double>::iterator it = elemental_error.begin(); it != elemental_error.end(); it++)
    {
      z2_global_sum += *it;
    }
    sprintf(filename,"%s/global_z2_error.dat", Doc_info.directory().c_str());
    some_file.open(filename);

    some_file << z2_global_sum << std::endl;
    some_file.close();

    oomph_info << "Global Z2 error: " << z2_global_sum << std::endl;
  }
  
  // increment the solution numbering
  Doc_info.number()++;
  
} // end_of_doc_solution

template<class ELEMENT>
void FibreInTwoDBoxProblem<ELEMENT>::impose_fake_singular_amplitude()
{
  // loop over the singular functions we're subtracting
  for(unsigned ising=0; ising<Nsingular_fct; ising++)
  {
    ScalableSingularityForNavierStokesElement<ELEMENT>* el_pt =
      dynamic_cast<ScalableSingularityForNavierStokesElement<ELEMENT>*>
      (Singular_fct_element_mesh_pt->element_pt(ising));
  
    // Change r_C so that C is assigned directly
    double imposed_amplitude = Global_Physical_Variables::singular_amplitude_for_debug;
    el_pt->impose_singular_fct_amplitude(imposed_amplitude);
  }
}

//== start of set_values_to_singular_solution ============================
/// Function to assign the singular solution to all nodes of the mesh
//========================================================================
template<class ELEMENT>
void FibreInTwoDBoxProblem<ELEMENT>::set_values_to_singular_solution()
{
  // get the number of nodes in the mesh
  unsigned nnode = Bulk_mesh_pt->nnode();
  
  for(unsigned i=0; i<nnode; i++)
  {
    // get a pointer to this node
    Node* node_pt = Bulk_mesh_pt->node_pt(i);

    // get the position of this node
    Vector<double> x(2, 0.0);
    x[0] = node_pt->x(0);
    x[1] = node_pt->x(1);

    // get the singular solution at this point
    Vector<double> u(3, 0.0);

    int boundary_id = -1;
    if(node_pt->is_on_boundary(Global_Physical_Variables::Inner_plate_boundary))
      boundary_id = Global_Physical_Variables::Inner_plate_boundary;
    else if(node_pt->is_on_boundary(Global_Physical_Variables::Inner_plate_boundary_upper))
      boundary_id = Global_Physical_Variables::Inner_plate_boundary_upper;

    // loop over the singular functions we're subtracting
    for(unsigned ising=0; ising<Nsingular_fct; ising++)
    {
      // call the singular elements analytic function
      Vector<double> u_ising =
	dynamic_cast<ScalableSingularityForNavierStokesElement<ELEMENT>*>(
	  Singular_fct_element_mesh_pt->element_pt(ising))->unscaled_singular_fct(x, boundary_id);

      // add this onto the total 
      for(unsigned j=0; j<Dim+1; j++)
      {
	u[j] += u_ising[j];
      }
    }
    
    // assign the velocities
    node_pt->set_value(0, u[0]);
    node_pt->set_value(1, u[1]);

    // this is a bit hacky but fuck it
    if(node_pt->nvalue() == 3 || node_pt->nvalue() == 5)
    {
      node_pt->set_value(2, u[2]);
    }
  }
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


//==start_of_main======================================================
/// Driver for backward step with impedance outflow bc
//=====================================================================
int main(int argc, char **argv)
{

  // Store command line arguments
  CommandLineArgs::setup(argc,argv);
  
  // results directory
  CommandLineArgs::specify_command_line_flag("--dir", &Global_Physical_Variables::dir);

  CommandLineArgs::specify_command_line_flag("-output_hi_res_soln");
  
  // specify non-default element area
  CommandLineArgs::specify_command_line_flag("--element_area", &Global_Physical_Variables::element_area);

  // specify the element area to use in the high-res regions around the plates edges
  CommandLineArgs::specify_command_line_flag(
    "--high_res_element_area", &Global_Physical_Variables::hi_res_element_area);
  
  // // specify the radius of the high resolution circular region around each edge of the plate
  // CommandLineArgs::specify_command_line_flag(
  //   "--high_res_region_radius", &Global_Physical_Variables::high_res_region_radius);
  
  // specify the x-component of the plate velocity
  CommandLineArgs::specify_command_line_flag("--velocity_x",
					     &Global_Physical_Variables::plate_velocity[0]);

  // specify the y-component of the plate velocity
  CommandLineArgs::specify_command_line_flag("--velocity_y",
					     &Global_Physical_Variables::plate_velocity[1]);

    // specify the radius of the plate curvature
  CommandLineArgs::specify_command_line_flag("--plate_radius_of_curvature",
					     &Global_Physical_Variables::plate_radius_of_curvature);
  
  // specify the the plate's angular velocity (right-handed, so positive is anticlockwise)
  CommandLineArgs::specify_command_line_flag("--angular_velocity_plate",
					     &Global_Physical_Variables::angular_velocity_plate);

    // specify the x-component of the velocity at the outer boundaries
  CommandLineArgs::specify_command_line_flag("--boundary_velocity_x", 
					     &Global_Physical_Variables::boundary_velocity[0]);
  
  // specify the y-component of the velocity at the outer boundaries
  CommandLineArgs::specify_command_line_flag("--boundary_velocity_y", 
					     &Global_Physical_Variables::boundary_velocity[1]);
  
  // specify the boundaries' angular velocity
  CommandLineArgs::specify_command_line_flag("--angular_velocity_boundary",
					     &Global_Physical_Variables::angular_velocity_boundary);
  
  // allow for no translational velocity of boundaries or plate, but rotation of just the boundary
  CommandLineArgs::specify_command_line_flag("--no_plate_rotation");
  
  // specify the x-coordinate of the centre of rotation of the plate
  CommandLineArgs::specify_command_line_flag("--centre_of_rotation_x",
					     &Global_Physical_Variables::centre_of_rotation[0]);
  
  // specify the y-coordinate of the centre of rotation of the plate
  CommandLineArgs::specify_command_line_flag("--centre_of_rotation_y",
					     &Global_Physical_Variables::centre_of_rotation[1]);

  // set the y-coordinate of the lower edge of the domain
  CommandLineArgs::specify_command_line_flag("--bottom_edge_y",
					     &Global_Physical_Variables::bottom_edge_y);

  // set the y-coordinate of the top edge of the domain
  CommandLineArgs::specify_command_line_flag("--top_edge_y",
					     &Global_Physical_Variables::top_edge_y);

  // set the x-coordinate of the left edge of the domain
  CommandLineArgs::specify_command_line_flag("--left_edge_x",
					     &Global_Physical_Variables::L_edge_x);

  // set the x-coordinate of the right edge of the domain
  CommandLineArgs::specify_command_line_flag("--right_edge_x",
					     &Global_Physical_Variables::R_edge_x);
  
  // just compute with pure finite elements rather than with the augmented singular elements
  CommandLineArgs::specify_command_line_flag("--dont_subtract_singularity");

  // assign the singular solution to all nodal values, for validation
  CommandLineArgs::specify_command_line_flag("--set_values_to_singular_solution");
  
  // output the Jacobian matrix in sparse-indexed form
  CommandLineArgs::specify_command_line_flag("--output_jacobian_sparse");

  // output the Z2 error at the locations of the singularities
  CommandLineArgs::specify_command_line_flag("--output_z2_error");
  
  // output description of dofs
  CommandLineArgs::specify_command_line_flag("--describe_dofs");

  // output the location of each node
  CommandLineArgs::specify_command_line_flag("--describe_nodes");
  
  // fake singular amplitude for debug
  CommandLineArgs::specify_command_line_flag("--set_sing_amplitude",
	&Global_Physical_Variables::singular_amplitude_for_debug);
  
  // use Lagrange multipliers to enforce dirichlet boundary conditions
  CommandLineArgs::specify_command_line_flag("--enforce_dirichlet_bcs_by_lagrange_multipliers");

  // switch off the hi-res regions for debug
  CommandLineArgs::specify_command_line_flag("--no_hi_res_regions");

  // switch from a flat plate to a circular plate
  CommandLineArgs::specify_command_line_flag("--circular_plate");
  
  // Parse command line
  CommandLineArgs::parse_and_assign(); 
  
  // Doc what has actually been specified on the command line
  CommandLineArgs::doc_specified_flags();

  // QUEHACERES for curved plates there seems to be a mismatch of 2e-14...
  ToleranceForVertexMismatchInPolygons::Tolerable_error = 3e-14;

  if (CommandLineArgs::command_line_flag_has_been_set("--no_hi_res_regions") )
  {
    Global_Physical_Variables::include_hi_res_regions = false;
    Global_Physical_Variables::Inner_plate_boundary_upper = 4;
  }
      
  // Build the problem with Triangular Taylor Hood elements
  FibreInTwoDBoxProblem<ProjectableTaylorHoodElement<MyTNavierStokesElement <2,3> > > problem;

  // set the global pointer
  Global_Physical_Variables::Plate_pt = problem.plate_pt();

  // are we imposing the amplitude directly and bypassing the calculation?
  if (CommandLineArgs::command_line_flag_has_been_set("--set_sing_amplitude") )
  {
    if (!CommandLineArgs::command_line_flag_has_been_set("--dont_subtract_singularity"))
    {
      problem.impose_fake_singular_amplitude();
    }
  }

  if (CommandLineArgs::command_line_flag_has_been_set("--set_values_to_singular_solution"))
  {
    problem.set_values_to_singular_solution();
  }

  // doc initial conditions
  problem.doc_solution();

  problem.max_residuals() = 1e3;
  
  // solve the bloody thing
  problem.newton_solve();

  // doc it
  problem.doc_solution();
  
  return 0;
}
