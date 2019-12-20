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
// Header file for elements that are used to ...hierher
#ifndef OOMPH_STOKES_SING_FACE_ELEMENTS_HEADER
#define OOMPH_STOKES_SING_FACE_ELEMENTS_HEADER

#include "navier_stokes.h"

namespace oomph
{

  //============================================================================
  // TemplateFreeScalableSingularityForNavierStokesElement defines the elements managing
  // the singular function : it is essentially a pointer to the singular function, 
  // its gradient and its amplitude
  //============================================================================
  class TemplateFreeScalableSingularityForNavierStokesElement :
    public virtual GeneralisedElement
  {
  public:

    typedef Vector<double>(*UnscaledSingSolnFctPt) (const Vector<double>& x, int boundary_id);


    typedef DenseMatrix<double>(*GradientOfUnscaledSingSolnFctPt) 
      (const Vector<double>& x, int boundary_id);
    
    ///Constructor
  TemplateFreeScalableSingularityForNavierStokesElement()    
    {
      //data to store amplitude
      add_internal_data(new Data(1));
    }

    ///Function to get pointer to unscaled version of singular function
    UnscaledSingSolnFctPt& unscaled_singular_fct_pt()
    {
      return Unscaled_singular_fct_pt;
    }

    ///Function to get pointer to unscaled version of gradient of singular function
    GradientOfUnscaledSingSolnFctPt& gradient_of_unscaled_singular_fct_pt() 
    {
      return Gradient_of_unscaled_singular_fct_pt;
    }

    ///Function to compute unscaled version of unscaled version
    Vector<double> unscaled_singular_fct(const Vector<double>& x, int boundary_id=-1) const
    {
      if(Unscaled_singular_fct_pt == 0)
      {
	return *(new Vector<double>(x.size(), 0.0));
      }
      return Unscaled_singular_fct_pt(x, boundary_id);
    }

    ///Compute unscaled version of gradient of singular function
    DenseMatrix<double> gradient_of_unscaled_singular_fct(const Vector<double>& x, int boundary_id=-1)
      const
    {
      DenseMatrix<double> grad;
      
      if(Gradient_of_unscaled_singular_fct_pt == 0)
      {
	return grad;
      }
      
      return Gradient_of_unscaled_singular_fct_pt(x, boundary_id);
    }

    ///Compute scaled version of singular function
    Vector<double> singular_fct(const Vector<double>& x, int boundary_id=-1) const
    {
      // get dimension of the problem; plus one because we want pressure as well
      // as the velocity components
      const unsigned Dim = x.size() + 1;

      // storage for the scaled basis functions
      Vector<double> scaled_singular_fct(Dim, 0.0);

      // get the unscaled functions
      Vector<double> unscaled_u_sing(Dim);
      unscaled_u_sing = unscaled_singular_fct(x, boundary_id);

      double amplitude = amplitude_of_singular_fct();
      
      // scale 'em
      for(unsigned i=0; i<Dim; i++)
      {
	scaled_singular_fct[i] = amplitude * unscaled_u_sing[i];
      }
      
      return scaled_singular_fct;
    }

    ///Compute scaled version of gradient of singular function
    DenseMatrix<double> gradient_of_singular_fct(const Vector<double>& x, int boundary_id=-1) const
    {
      DenseMatrix<double> grad(gradient_of_unscaled_singular_fct(x, boundary_id));
      
      const unsigned n = grad.nrow();
      const unsigned m = grad.ncol();
      
      for(unsigned i=0; i<n; i++)
      {
	for(unsigned j=0; j<m; j++)
	{
	  grad(i,j) *= amplitude_of_singular_fct();
	}
      }
      return grad;
    }

    ///Access the amplitude of the singular function
    double amplitude_of_singular_fct() const
    {
      return data_that_stores_amplitude_of_singular_fct()
        ->value(index_of_value_that_stores_amplitude_of_singular_fct());
    }

    ///Set the amplitude of thz singular function
    void set_amplitude_of_singular_fct(const double& value)
    {
      data_that_stores_amplitude_of_singular_fct()
        ->set_value(index_of_value_that_stores_amplitude_of_singular_fct(),value);
    }

    ///pin amplitude of singular function
    void pin_amplitude_of_singular_fct()
    {
      data_that_stores_amplitude_of_singular_fct()
        ->pin(index_of_value_that_stores_amplitude_of_singular_fct());
    }

    ///Pointer to data that stores the amplitude of singular function
    Data* data_that_stores_amplitude_of_singular_fct() const
    {
      return internal_data_pt(0);
    }

    ///Gives the index of the amplitude value : default is 0
    unsigned index_of_value_that_stores_amplitude_of_singular_fct() const 
    {
      return 0;
    }
    
  private:

    ///Pointer to singular function
    UnscaledSingSolnFctPt Unscaled_singular_fct_pt;

    ///Pointer to gradient of singular funcion
    GradientOfUnscaledSingSolnFctPt Gradient_of_unscaled_singular_fct_pt;
  };

  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////

  
  //===========================================================================
  /// NavierStokesWithSingularityBoundaryIntegralFaceElement is a class of face elements 
  ///used to compute the contribution to the residuals from the the singular function
  //===========================================================================
  template <class ELEMENT>
    class NavierStokesWithSingularityBoundaryIntegralFaceElement :
    public virtual FaceGeometry<ELEMENT>, public virtual FaceElement
  {
      
  public:

    /// \short Function pointer to the "exact" non-singular function fct(x,u,grad u)
    typedef void (*ExactNonSingularFctPt)(const Vector<double>& x, Vector<double>& u,
					  DenseMatrix<double>& grad_u);
      
    /// \short Constructor, takes the pointer to the "bulk" element and the 
    /// index of the face to which the element is attached.
    NavierStokesWithSingularityBoundaryIntegralFaceElement(FiniteElement* const& bulk_el_pt, 
							   const int& face_index,
							   const int& boundary_id=0);

    ///\short  Broken empty constructor
    NavierStokesWithSingularityBoundaryIntegralFaceElement()
    {
      throw OomphLibError(
	"Don't call empty constructor for NavierStokesWithSingularityBoundaryIntegralFaceElement",
	OOMPH_CURRENT_FUNCTION,
	OOMPH_EXCEPTION_LOCATION);
    }

    /// Broken copy constructor
    NavierStokesWithSingularityBoundaryIntegralFaceElement(
      const NavierStokesWithSingularityBoundaryIntegralFaceElement& dummy) 
    { 
      BrokenCopy::broken_copy("NavierStokesWithSingularityBoundaryIntegralFaceElement");
    } 
      
    /// Broken assignment operator
    void operator=(const NavierStokesWithSingularityBoundaryIntegralFaceElement&) 
      {
	BrokenCopy::broken_assign("NavierStokesWithSingularityBoundaryIntegralFaceElement");
      }

    /// \short Specify the value of nodal zeta from the face geometry
    /// The "global" intrinsic coordinate of the element when
    /// viewed as part of a geometric object should be given by
    /// the FaceElement representation, by default (needed to break
    /// indeterminacy if bulk element is SolidElement)
    double zeta_nodal(const unsigned &n, const unsigned &k,           
		      const unsigned &i) const 
    {
      return FaceElement::zeta_nodal(n,k,i);
    }

    /// Pointer to element that computes singular function related stuff
    TemplateFreeScalableSingularityForNavierStokesElement*& navier_stokes_sing_el_pt()
    { 
      return Navier_stokes_sing_el_pt;
    } 
 
    /// Add the element's contribution to its residual vector
    inline void fill_in_contribution_to_residuals(Vector<double> &residuals)
    {
      //Call the generic residuals function with flag set to 0
      //using a dummy matrix argument
      fill_in_generic_residual_contribution_to_singular_amplitude(
	residuals,GeneralisedElement::Dummy_matrix,0);
    }

#ifndef USE_FD_JACOBIAN
    /// \short Add the element's contribution to its residual vector and its
    /// Jacobian matrix
    inline void fill_in_contribution_to_jacobian(Vector<double> &residuals,
						 DenseMatrix<double> &jacobian)
    {
      //Call the generic routine with the flag set to 1
      fill_in_generic_residual_contribution_to_singular_amplitude(residuals,jacobian,1);
    }
#endif
    
    /// Output function -- forward to broken version in FiniteElement
    /// until somebody decides what exactly they want to plot here...
    void output(std::ostream &outfile)
    {
      FiniteElement::output(outfile);
    }

    /// Output with various contributions
    void output(std::ostream &outfile, const unsigned &nplot)
    {
      //Vector of local coordinates
      Vector<double> s(Dim-1);
   
      // Tecplot header info
      outfile << this->tecplot_zone_string(nplot);
   
      // Loop over plot points
      unsigned num_plot_points = this->nplot_points(nplot);
      for (unsigned iplot=0; iplot<num_plot_points; iplot++)
      {
	// Get local coordinates of plot point
	this->get_s_plot(iplot,nplot,s);
     
	Vector<double> x(Dim);
	for(unsigned i=0; i<Dim; i++) 
	{
	  x[i]=this->interpolated_x(s,i);
	  outfile << x[i] << " ";
	}

	// Compute outer unit normal at the specified local coordinate
	Vector<double> unit_normal(Dim);
	outer_unit_normal(s,unit_normal);

	outfile << unit_normal[0] << " " << unit_normal[1] << " " << endl;

	// QUEHACERES output something useful here
      }
   
      // Write tecplot footer (e.g. FE connectivity lists)
      this->write_tecplot_zone_footer(outfile,nplot);
   
    }


    /// C-style output function -- forward to broken version in FiniteElement
    /// until somebody decides what exactly they want to plot here...
    void output(FILE* file_pt)
    {
      FiniteElement::output(file_pt);
    }

    /// \short C-style output function -- forward to broken version in 
    /// FiniteElement until somebody decides what exactly they want to plot 
    /// here...
    void output(FILE* file_pt, const unsigned &n_plot)
    {
      FiniteElement::output(file_pt,n_plot);
    }

    /// \short Compute this element's contribution to the integral that determines C
    double get_contribution_integral();
 
    /// Pointer to exact non singular fct (and its gradient) only used
    /// to validate the computation of the integral. Ignored if null 
    /// which is the default
    // hierher make read only and use set/unset fct to enable/disable;
    // currently we'd have to reset this to null!
    ExactNonSingularFctPt& exact_non_singular_fct_pt()
    {
      return Exact_non_singular_fct_pt;
    }

    // get the ID of the boundary this face element sits on
    int boundary_id()
    {
      return Boundary_id;
    }
      
  protected:

    /// \short Function to compute the shape and test functions and to return 
    /// the Jacobian of mapping between local and global (Eulerian)
    /// coordinates
    inline double shape_and_test(const Vector<double>& s, Shape& psi, Shape& test)
      const
    {
      //Find number of nodes
      unsigned n_node = nnode();

      //Get the shape functions
      shape(s,psi);

      //Set the test functions to be the same as the shape functions
      for(unsigned i=0; i<n_node; i++)
      {
	test[i] = psi[i];
      }

      //Return the value of the jacobian
      return J_eulerian(s);
    }

    /// \short Function to compute the shape and test functions and to return 
    /// the Jacobian of mapping between local and global (Eulerian)
    /// coordinates
    inline double shape_and_test_at_knot(const unsigned& ipt,
					 Shape& psi, Shape& test)
      const
    {
      //Find number of nodes
      unsigned n_node = nnode();

      //Get the shape functions
      shape_at_knot(ipt,psi);

      //Set the test functions to be the same as the shape functions
      for(unsigned i=0; i<n_node; i++)
      {
	test[i] = psi[i];
      }

      //Return the value of the jacobian
      return J_eulerian_at_knot(ipt);
    }

  private:

    // ID of the boundary that this face element sits on
    int Boundary_id;
      
    /// \short Add the element's contribution to its residual vector.
    /// flag=1(or 0): do (or don't) compute the contribution to the
    /// Jacobian as well. 
    void fill_in_generic_residual_contribution_to_singular_amplitude(
      Vector<double> &residuals, DenseMatrix<double> &jacobian, 
      const unsigned& flag);
 
    ///The spatial dimension of the problem
    unsigned Dim;

    ///The index at which the FE pressure is stored at the nodes
    unsigned P_index_nst;

    /// Pointer to element that stores the singular fcts etc.
    TemplateFreeScalableSingularityForNavierStokesElement* Navier_stokes_sing_el_pt;

    /// Pointer to exact non singular fct (and its gradient) only used
    /// to validate the computation of the integral. Ignored if null 
    /// which is the default
    ExactNonSingularFctPt Exact_non_singular_fct_pt;
  }; 

  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  //===========================================================================
  /// Constructor, takes the pointer to the "bulk" element, the 
  /// index of the fixed local coordinate and its value represented
  /// by an integer (+/- 1), indicating that the face is located
  /// at the max. or min. value of the "fixed" local coordinate
  /// in the bulk element.
  //===========================================================================
  template<class ELEMENT>
    NavierStokesWithSingularityBoundaryIntegralFaceElement<ELEMENT>::
    NavierStokesWithSingularityBoundaryIntegralFaceElement(FiniteElement* const& bulk_el_pt, 
							   const int& face_index,
							   const int& boundary_id) : 
  FaceGeometry<ELEMENT>(), FaceElement()
  { 

    // Null out the fct pointer so integral is computed with
    // actual finite element representation of the non singular
    // part
    Exact_non_singular_fct_pt = 0;

    // set the ID of the boundary that this face element sits on
    Boundary_id = boundary_id;
    
    // Let the bulk element build the FaceElement, i.e. setup the pointers 
    // to its nodes (by referring to the appropriate nodes in the bulk
    // element), etc.
    bulk_el_pt->build_face_element(face_index,this);

#ifdef PARANOID
    {
      //Check that the element is not a refineable 3d element
      ELEMENT* elem_pt = dynamic_cast<ELEMENT*>(bulk_el_pt);
      //If it's three-d
      if(elem_pt->dim()==3)
      {
	//Is it refineable
	RefineableElement* ref_el_pt=dynamic_cast<RefineableElement*>(elem_pt);
	if(ref_el_pt!=0)
	{
	  if (this->has_hanging_nodes())
	  {
	    throw OomphLibError(
	      "This traction element will not work correctly if nodes are hanging\n",
	      OOMPH_CURRENT_FUNCTION,
	      OOMPH_EXCEPTION_LOCATION);
	  }
	}
      }
    }
#endif

    // Initialising the pointer to the singularity function
    this->Navier_stokes_sing_el_pt = 0;
 
    // Extract the dimension of the problem from the dimension of 
    // the first node
    Dim = this->node_pt(0)->ndim();

    // Set up P_index_nst. Initialise to zero, which probably won't change
    // in most cases, oh well, the price we pay for generality
    P_index_nst = 0;

    // Cast to the appropriate NavierStokesEquation so that we can
    // find the index at which the variable is stored
    // We assume that the dimension of the full problem is the same
    // as the dimension of the node, if this is not the case you will have
    // to write custom elements, sorry
    switch(Dim)
    {
      //One dimensional problem
      case 1:
      {
	NavierStokesEquations<1>* eqn_pt = 
	  dynamic_cast<NavierStokesEquations<1>*>(bulk_el_pt);
	//If the cast has failed die
	if(eqn_pt == 0)
	{
	  std::string error_string =
	    "Bulk element must inherit from NavierStokesEquations.";
	  error_string += 
	    "Nodes are one dimensional, but cannot cast the bulk element to\n";
	  error_string += "NavierStokesEquations<1>\n.";
	  error_string += 
	    "If you desire this functionality, you must implement it yourself\n";
       
	  throw OomphLibError(error_string,
			      OOMPH_CURRENT_FUNCTION,
			      OOMPH_EXCEPTION_LOCATION);
	}
	//Otherwise read out the value
	else
	{
	  //Read the index from the (cast) bulk element
	  P_index_nst = eqn_pt->p_nodal_index_nst();
	}
      }
      break;
    
      //Two dimensional problem
      case 2:
      {
	NavierStokesEquations<2>* eqn_pt = 
	  dynamic_cast<NavierStokesEquations<2>*>(bulk_el_pt);
	//If the cast has failed die
	if(eqn_pt==0)
	{
	  std::string error_string =
	    "Bulk element must inherit from NavierStokesEquations.";
	  error_string += 
	    "Nodes are two dimensional, but cannot cast the bulk element to\n";
	  error_string += "NavierStokesEquations<2>\n.";
	  error_string += 
	    "If you desire this functionality, you must implement it yourself\n";
       
	  throw OomphLibError(error_string,
			      OOMPH_CURRENT_FUNCTION,
			      OOMPH_EXCEPTION_LOCATION);
	}
	else
	{
	  //Read the index from the (cast) bulk element.
	  P_index_nst = eqn_pt->p_nodal_index_nst();
	}
      }
      break;
    
      //Three dimensional problem
      case 3:
      {
	NavierStokesEquations<3>* eqn_pt = 
	  dynamic_cast<NavierStokesEquations<3>*>(bulk_el_pt);
	//If the cast has failed die
	if(eqn_pt==0)
	{
	  std::string error_string =
	    "Bulk element must inherit from NavierStokesEquations.";
	  error_string += 
	    "Nodes are three dimensional, but cannot cast the bulk element to\n";
	  error_string += "NavierStokesEquations<3>\n.";
	  error_string += 
	    "If you desire this functionality, you must implement it yourself\n";
       
	  throw OomphLibError(error_string,
			      OOMPH_CURRENT_FUNCTION,
			      OOMPH_EXCEPTION_LOCATION);
       
	}
	else
	{
	  //Read the index from the (cast) bulk element.
	  P_index_nst = eqn_pt->p_nodal_index_nst();
	}
      }
      break;

      //Any other case is an error
      default:
	std::ostringstream error_stream; 
	error_stream <<  "Dimension of node is " << Dim 
		     << ". It should be 1,2, or 3!" << std::endl;
     
	throw OomphLibError(error_stream.str(),
			    OOMPH_CURRENT_FUNCTION,
			    OOMPH_EXCEPTION_LOCATION);
	break;
    }
  }

  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  //===========================================================================
  /// Compute the element's residual vector and the (zero) Jacobian matrix.
  //===========================================================================
  template<class ELEMENT>
    void NavierStokesWithSingularityBoundaryIntegralFaceElement<ELEMENT>::
    fill_in_generic_residual_contribution_to_singular_amplitude(
      Vector<double> &residuals, DenseMatrix<double> &jacobian, 
      const unsigned& flag)
  {

    // hierher populate when contribution is split up
    oomph_info << "This shouldn't be called at the moment\n";
    abort();
  }

  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  
  //===========================================================================
  /// Calculate the contribution of the face element to the integral that
  /// determines the amplitude, via the Lorentz reciprocity theorem. 
  //===========================================================================
  template<class ELEMENT>
    double NavierStokesWithSingularityBoundaryIntegralFaceElement<ELEMENT>::
    get_contribution_integral()
  {
    //Find out how many nodes there are
    const unsigned n_node = nnode();
  
    //Set up memory for the shape and test functions
    Shape psif(n_node), testf(n_node);
 
    //Set the value of Nintpt
    const unsigned n_intpt = integral_pt()->nweight();
 
    // Set the Vector to hold local coordinates (in this face element, not the
    // bulk element this is attached to)
    Vector<double> s(Dim-1);
 
    // Saves result of integration
    double integral_result = 0.0;

    //Loop over the integration points
    //--------------------------------
    for(unsigned ipt=0; ipt<n_intpt; ipt++)
    {

      //Assign values of s
      for(unsigned i=0; i<(Dim-1); i++)
      {
	s[i] = integral_pt()->knot(ipt,i);
      }
   
      //Get the integral weight
      double w = integral_pt()->weight(ipt);
   
      //Find the shape and test functions and return the Jacobian
      //of the mapping
      double J = shape_and_test(s, psif, testf);
   
      //Premultiply the weights and the Jacobian
      double W = w*J;

      // compute outer normal unit vector
      Vector<double> unit_normal(Dim);
      outer_unit_normal(s, unit_normal);

      // local coordinates in bulk element this face element is attached to
      Vector<double> s_bulk(Dim);

      // global coordinates
      Vector<double> x(Dim); 

      // get global coordinates
      for(unsigned i=0; i<Dim; i++)
      { 
	x[i] = this->interpolated_x(s,i); 
      } 
   
      // Get gradient of scaled/unscaled singular velocity functions
      DenseMatrix<double> dudx_sing(Dim, Dim);
      DenseMatrix<double> dudx_sing_unscaled(Dim, Dim);
      dudx_sing          = Navier_stokes_sing_el_pt->gradient_of_singular_fct(x);
      dudx_sing_unscaled = Navier_stokes_sing_el_pt->gradient_of_unscaled_singular_fct(x);
      
      // Get the values of the singular functions at our current location
      Vector<double> u_sing(Dim+1);
      Vector<double> u_sing_unscaled(Dim+1);
      u_sing          = Navier_stokes_sing_el_pt->singular_fct(x);
      u_sing_unscaled = Navier_stokes_sing_el_pt->unscaled_singular_fct(x);

      // get singular pressure
      double p_sing          = u_sing[P_index_nst];
      double p_sing_unscaled = u_sing_unscaled[P_index_nst];
      
      // shorthand
      ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(this->bulk_element_pt());

      // compute the singular contribution to the strain-rate
      DenseMatrix<double>strain_rate_sing(Dim, Dim, 0.0);
      DenseMatrix<double>strain_rate_sing_unscaled(Dim, Dim, 0.0);
      
      for (unsigned i=0; i<Dim; i++)
      {
	for(unsigned j=0; j<Dim; j++)
	{
	  strain_rate_sing(i,j) = 0.5*(dudx_sing(i,j) + dudx_sing(j,i));
	  
	  strain_rate_sing_unscaled(i,j) =
	    0.5*(dudx_sing_unscaled(i,j) + dudx_sing_unscaled(j,i));
	}
      }
	
      // get contribution of singular pressure and singular velocity gradients to stress tensor
      DenseMatrix<double> stress_sing(Dim, Dim);
      DenseMatrix<double> stress_sing_unscaled(Dim, Dim);
      stress_sing = (*bulk_el_pt->stress_fct_pt())(strain_rate_sing, p_sing);
      stress_sing_unscaled = (*bulk_el_pt->stress_fct_pt())(strain_rate_sing_unscaled, p_sing_unscaled);
      
      // Get the local bulk coordinates    
      s_bulk = local_coordinate_in_bulk(s);
      
      Vector<double> u_fe(Dim);
      
      // Get FE part of the velocity
      for(unsigned i=0; i<Dim; i++)
      {
	u_fe[i] = bulk_el_pt->interpolated_u_nst(s_bulk, i);
      }
      
      // get FE part of pressure
      double p_fe = bulk_el_pt->interpolated_p_nst(s_bulk);

      // get FE part of velocity gradient tensor
      DenseMatrix<double> dudx_fe(Dim, Dim);
      
      for(unsigned i=0; i<Dim; i++)
      {
	for(unsigned j=0; j<Dim; j++)
	{
	  // get derivative du_i/dx_j
	  dudx_fe(i,j) = bulk_el_pt->interpolated_dudx_nst(s_bulk, i, j);
	}
      }
      
      // get the FE strain rate 1/2(du_i/dx_j + du_j/dx_i)
      DenseMatrix<double> strain_rate_fe(Dim, Dim);
      
      bulk_el_pt->strain_rate(s_bulk, strain_rate_fe);
      
      // FE part of the stress
      DenseMatrix<double> stress_fe(Dim, Dim);

      // compute it from consitutive equation
      stress_fe = (*bulk_el_pt->stress_fct_pt())(strain_rate_fe, p_fe);
            
      // get FE part of the traction
      Vector<double> traction_fe(Dim);

      // get FE traction from bulk element
      bulk_el_pt->get_traction(s_bulk, unit_normal, traction_fe);

      // ================================================
      // Now we compute the contibution of the integral

      for(unsigned i=0; i<Dim; i++)
      {
	for(unsigned j=0; j<Dim; j++)
	{
	  // Lorentz reciprocity theorem
	  integral_result += W * unit_normal[j] * (stress_fe(i,j) * u_sing_unscaled[i]
						   - stress_sing_unscaled(i,j) * u_fe[i] );
	}
      }      
    } // end loop over integration points
        
    return integral_result;
  }




  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////


  //======================================================================
  /// \short Class for elements that handle singularities
  /// in the Navier-Stokes equations. Templated by bulk element within
  /// which we impose regularity on the FE solution by insisting that
  /// the slope of the solution at a specified local coordinate, and in
  /// in a specified direction is zero. Nodal values of that element
  /// become external data for the current element whose equation
  /// (zero slope of the FE solution, as discussed) determines the 
  /// amplitude of the singular function.
  //======================================================================
  template<class BULK_ELEMENT> 
    class ScalableSingularityForNavierStokesElement : 
    public virtual TemplateFreeScalableSingularityForNavierStokesElement
  {
   
  public:
   
    /// Constructor
  ScalableSingularityForNavierStokesElement() :
    Bulk_element_pt(0), Face_element_mesh_pt(0), 
      Impose_singular_fct_amplitude(false)
      { }
    
    /// Set pointer to mesh containing the FaceElements (and flush
    /// the previous ones first!)
    void set_mesh_of_face_elements(Mesh* const& face_mesh_pt)
    {
      Face_element_mesh_pt = face_mesh_pt;
      flush_external_data();
      
      unsigned nel = face_mesh_pt->nelement();
      
      oomph_info << "number of face elements used to compute C: "
		 << nel << std::endl;
      
      for (unsigned e=0; e<nel; e++)
      {	
	FiniteElement* el_pt =
	  dynamic_cast<NavierStokesWithSingularityBoundaryIntegralFaceElement<BULK_ELEMENT>*>(
	    face_mesh_pt->element_pt(e))->bulk_element_pt();
	
	unsigned nnod = el_pt->nnode();
	
	for (unsigned j=0; j<nnod; j++)
	{
	  add_external_data(el_pt->node_pt(j));
	}
      }
    }

    /// Call this to bypass the correct computation of the
    /// residual and replace it by r_c = C-ampl
    void impose_singular_fct_amplitude(double const& ampl)
    {
      Impose_singular_fct_amplitude = true;
      Imposed_amplitude = ampl;
    } 

    /// Reset to compute r_c properly via integral
    void dont_impose_singular_fct_amplitude()
    {
      Impose_singular_fct_amplitude = false;
    } 

    /// Add the element's contribution to its residual vector
    inline void fill_in_contribution_to_residuals(Vector<double> &residuals)
    {
      //Call the generic residuals function with flag set to 0
      //using a dummy matrix argument
      fill_in_generic_residual_contribution_navier_stokes_sing_fct(
	residuals,GeneralisedElement::Dummy_matrix, 0);
    }

    // function to return whether we're imposing the amplitude or not
    bool is_singular_fct_amplitude_imposed()
    {
      return Impose_singular_fct_amplitude;
    }
      
  private:

    /// Add the element's contribution to its residual vector
    inline void fill_in_generic_residual_contribution_navier_stokes_sing_fct(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      const unsigned& flag)
    {
      
      if (ndof() == 0)
      {
	return;
      }

#ifdef PARANOID
      // hierher paranoid check null pointers and zero sized vectors    
#endif
      
      // Get eqn number of residual that determines C
      int local_eqn_c = internal_local_eqn(0, 0);
      
      if (local_eqn_c >= 0)
      {
	// Bypass actual computation?
	if (Impose_singular_fct_amplitude)
	{
	  residuals[local_eqn_c] = this->amplitude_of_singular_fct() - Imposed_amplitude;
	}	
	// Do it properly
	else
	{
	  unsigned n_element = Face_element_mesh_pt->nelement();
	  for(unsigned e = 0; e<n_element; e++)
	  {
	    // get a pointer to this boundary face element
	    NavierStokesWithSingularityBoundaryIntegralFaceElement
	      <BULK_ELEMENT>* face_elem_pt =
	      dynamic_cast<NavierStokesWithSingularityBoundaryIntegralFaceElement
	      <BULK_ELEMENT>*>( Face_element_mesh_pt->finite_element_pt(e) );

	    // add it to the r_c residual equation
	    residuals[local_eqn_c] += face_elem_pt->get_contribution_integral();
	  }	  
	}
      }
    }

  private:  
  
    /// Pointer to bulk element where FE solution is regularised
    BULK_ELEMENT* Bulk_element_pt;

    /// Pointer to mesh of face elements that contribute to the surface
    /// integral that determines the amplitude of the unkown function
    Mesh* Face_element_mesh_pt;
  
    /// Imposed amplitude (only used if Impose_singular_fct_amplitude=true)
    double Imposed_amplitude;  

    /// \short Boolean to bypass the correct computation of the
    /// residual and replace it by r_c = C-ampl
    bool Impose_singular_fct_amplitude;

    // QUEHACERES for debug, contribution to the r_c integral from each boundary
    Vector<double> Integral_contribution_from_boundary;
  };




  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  
  //======================================================================
  /// \short A class for elements that imposes Dirichlet boundary 
  /// conditions on complete solution (such that u_fe + C u_sing = u_bc) using a
  /// Lagrange multiplier. Thus the element introduce an additional
  /// unknown at the nodes it's attached to. C and u_sing are specified
  /// via a ScalableSingularityForNavierStokesElement.
  //======================================================================
  template <class ELEMENT>
    class NavierStokesWithSingularityBCFaceElement : 
    public virtual FaceGeometry<ELEMENT>, 
    public virtual FaceElement 
    {
 
    public:

      /// \short Constructor, takes the pointer to the "bulk" element and the 
      /// index of the face to which the element is attached. Optional final
      /// arg is the identifier for the additional unknowns multiplier
      NavierStokesWithSingularityBCFaceElement(FiniteElement* const &bulk_el_pt, 
					       const int& face_index,
					       const unsigned &id = 0); 
  
      ///\short  Broken empty constructor
      NavierStokesWithSingularityBCFaceElement()
      {
	throw OomphLibError(
	  "Don't call empty constructor for NavierStokesWithSingularityBCFaceElement",
	  OOMPH_CURRENT_FUNCTION,
	  OOMPH_EXCEPTION_LOCATION);
      }
  
      /// Broken copy constructor
      NavierStokesWithSingularityBCFaceElement(
	const NavierStokesWithSingularityBCFaceElement& dummy) 
      { 
	BrokenCopy::broken_copy("NavierStokesWithSingularityBCFaceElement");
      } 
  
      /// Broken assignment operator
      void operator=(const NavierStokesWithSingularityBCFaceElement&) 
	{
	  BrokenCopy::broken_assign("NavierStokesWithSingularityBCFaceElement");
	}
  
      /// \short Specify the value of nodal zeta from the face geometry
      /// The "global" intrinsic coordinate of the element when
      /// viewed as part of a geometric object should be given by
      /// the FaceElement representation, by default (needed to break
      /// indeterminacy if bulk element is SolidElement)
      double zeta_nodal(const unsigned &n, const unsigned &k,           
			const unsigned &i) const 
      {
	return FaceElement::zeta_nodal(n,k,i);
      }

      /// Pointer to element that handles the ith singular fct
      ScalableSingularityForNavierStokesElement<ELEMENT>*
	navier_stokes_sing_el_pt(const unsigned& i) const
      {
	return Navier_stokes_sing_el_pt[i];
      }

      /// \short Set pointer to element that stores singular fct. Data that stores
      /// the amplitude of the singular fct and its index is retrieved from
      /// that element so the Data can be used as external Data in this
      /// element.
      void set_navier_stokes_sing_el_pt(
	Vector<ScalableSingularityForNavierStokesElement<ELEMENT>*>navier_stokes_sing_el_pt) 
      {
	// set the number of singular functions
	Nsingular_fct = navier_stokes_sing_el_pt.size();

	// make sure we've got enough space
	Navier_stokes_sing_el_pt.resize(Nsingular_fct);
	C_external_data_index.resize(Nsingular_fct);
	C_external_data_value_index.resize(Nsingular_fct);

	// loop over the singular functions and add their amplitudes as external data
	for(unsigned ising=0; ising<Nsingular_fct; ising++)
	{
	  Navier_stokes_sing_el_pt[ising] = navier_stokes_sing_el_pt[ising];
	  C_external_data_index[ising] = add_external_data(
	    navier_stokes_sing_el_pt[ising]->data_that_stores_amplitude_of_singular_fct());
	  C_external_data_value_index[ising] =
	    navier_stokes_sing_el_pt[ising]->index_of_value_that_stores_amplitude_of_singular_fct();
	}
      }


      /// Add the element's contribution to its residual vector
      inline void fill_in_contribution_to_residuals(Vector<double> &residuals)
      {
	//Call the generic residuals function with flag set to 0
	//using a dummy matrix argument
	fill_in_generic_residual_contribution_navier_stokes_sing(
	  residuals, GeneralisedElement::Dummy_matrix, 0);
      }

#ifndef USE_FD_JACOBIAN
      /// \short Add the element's contribution to its residual vector and its
      /// Jacobian matrix
      inline void fill_in_contribution_to_jacobian(Vector<double> &residuals,
      						   DenseMatrix<double> &jacobian)
      {
      	//Call the generic routine with the flag set to 1
      	fill_in_generic_residual_contribution_navier_stokes_sing(residuals, jacobian, 1);
      }
#endif

      /// Output function
      void output(std::ostream &outfile)
      {
	const unsigned n_plot=5;
	output(outfile,n_plot);
      }

      /// \short Output function
      void output(std::ostream &outfile, const unsigned &nplot)
      {
	//oomph_info << "hierher need to update output fct" << std::endl;
	//Vector of local coordinates
	Vector<double> s(Dim-1);
   
	// Tecplot header info
	outfile << this->tecplot_zone_string(nplot);
   
	// Loop over plot points
	unsigned num_plot_points=this->nplot_points(nplot);
	for (unsigned iplot=0;iplot<num_plot_points;iplot++)
	{
	  // Get local coordinates of plot point
	  this->get_s_plot(iplot,nplot,s);
     
	  Vector<double> x(Dim);
	  for(unsigned i=0; i<Dim; i++) 
	  {
	    x[i]=this->interpolated_x(s,i);
	    outfile << x[i] << " ";
	  }
	  outfile << endl;
	}
	return;
      }

      /// \short Provide nodal values of desired boundary values.
      /// They're imposed by Lagrange multipliers.
      void set_nodal_boundary_values(const DenseMatrix<double>& nodal_boundary_value)
      {
#ifdef PARANOID
	if (nodal_boundary_value.nrow() != nnode())
	{
	  std::stringstream error;
	  error << "nodel_boundary_value is a matrix with " 
		<< nodal_boundary_value.nrow() 
		<< " rows, but should have the same number of rows as the number of nodes, "
		<< nnode();
	  throw OomphLibError(error.str(),
			      OOMPH_CURRENT_FUNCTION,
			      OOMPH_EXCEPTION_LOCATION);
	}
#endif
	Nodal_boundary_value = nodal_boundary_value;
      }

      /// Pin Lagrange multiplier associated with ith coordinate at specified local node
      void pin_lagrange_multiplier_at_specified_local_node(const unsigned& j,
							   const unsigned& i,
							   const int& id = -1)
      {
	// get the face IDs map for this node
	map<unsigned, unsigned> map_l = *(
	  dynamic_cast<BoundaryNodeBase*>(node_pt(j))->
	  index_of_first_value_assigned_by_face_element_pt() );

	unsigned lambda_index;

	// if no id specified, just take the index for the first (and probably only)
	// boundary in the map
	if(id == -1)
	{
	  lambda_index = map_l.begin()->second;
	}	
	else
	{
	  // otherwise, get the nodal index for the specified boundary ID
	  lambda_index = map_l[id];
	}
	node_pt(j)->pin(lambda_index+i);
      }

      /// Unpin ith component of FE part of the solution at specified local node
      void unpin_u_fe_at_specified_local_node(const unsigned& j, const unsigned& i)
      {   
	node_pt(j)->unpin(i);	  	
      }

      /// C-style output function -- forward to broken version in FiniteElement
      /// until somebody decides what exactly they want to plot here...
      void output(FILE* file_pt)
      {
	FiniteElement::output(file_pt);
      }

      /// \short C-style output function -- forward to broken version in 
      /// FiniteElement until somebody decides what exactly they want to plot 
      /// here...
      void output(FILE* file_pt, const unsigned& n_plot)
      {
	FiniteElement::output(file_pt, n_plot);
      }

      // QUEHACERES for debug, output the value of the Lagrange multipliers on the Dirichlet
      // boundaries and the associated stress
      void output_lagrange_multiplier_and_stress(std::ostream& outfile,
						 const unsigned& nplot )
      {
	// pointer to the bulk element this face element is attached to
	ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(this->bulk_element_pt());
      
	// get the number of nodes in this face element
	unsigned nnode = this->nnode();
    
	//Set up memory for the shape and test functions
	Shape psi(nnode), test(nnode);
    
	unsigned num_plot_points = this->nplot_points(nplot);
	for (unsigned iplot=0; iplot < num_plot_points; iplot++)
	{
	  // Get local coordinates of this plot point
	  Vector<double> s(this->dim());      
	  this->get_s_plot(iplot, nplot, s);
      
	  //Find the shape and test functions and return the Jacobian
	  //of the mapping
	  double J = this->shape_and_test(s, psi, test);
         
	  // local coordinates in bulk element this face element is attached to
	  Vector<double> s_bulk(Dim);

	  // Get the local bulk coordinates    
	  s_bulk = local_coordinate_in_bulk(s);

	  Vector<double> x(Dim);

	  for(unsigned i=0; i<Dim; i++)
	  {
	    x[i] = this->interpolated_x(s,i);
	  }
	  
	  // Lagrange multipliers
	  Vector<double> lambda(Dim);
      
	  unsigned nnode = this->nnode();
	  for(unsigned j=0; j<nnode; j++)
	  {
	    // get the map which gives the starting nodal index for
	    // the Lagrange multipliers associated with each boundary ID
	    std::map<unsigned, unsigned> first_index = *(
	      dynamic_cast<BoundaryNodeBase*>(this->node_pt(j))->
	      index_of_first_value_assigned_by_face_element_pt() );
	
	    for(unsigned i=0; i<Dim; i++)
	    {
	      unsigned lambda_index = first_index[Boundary_id] + i;
	  
	      // get the Lagrange multiplier
	      lambda[i] += this->nodal_value(j, lambda_index) * psi[j];
	    }
	  }

	  // get FE part of pressure
	  double p_fe = bulk_el_pt->interpolated_p_nst(s_bulk);

	  // total pressure
	  double p_total = p_fe;

	  // get the singular parts of the pressure
	  Vector<double> p_sing(Nsingular_fct, 0.0);
	  
	  for(unsigned ising=0; ising<Nsingular_fct; ising++)
	  {
	    // Get the values of the scaled singular functions at our current location
	    Vector<double> u_sing(Dim+1);
	    u_sing = Navier_stokes_sing_el_pt[ising]->singular_fct(x);
	    
	    // get scaled singular pressure
	    p_sing[ising] = u_sing[P_index_nst];

	    // add to the total pressure
	    p_total += p_sing[ising];
	  }
	  
	  // get the FE strain rate 1/2(du_i/dx_j + du_j/dx_i)
	  DenseMatrix<double> strain_rate_fe(Dim, Dim);      
	  bulk_el_pt->strain_rate(s_bulk, strain_rate_fe);

	  // get the stress from the constitutive law
	  DenseMatrix<double> stress_fe(Dim, Dim);
	  stress_fe = (*bulk_el_pt->stress_fct_pt())(strain_rate_fe, p_fe);
	  	  

	  // ====================================================
	  // output stuff
	  // ====================================================
	  
	  // output the coordinates of this point
	  for(unsigned i=0; i<Dim; i++)
	  {	    
	    outfile << x[i] << " ";
	  }
	  
	  // output the Lagrange multipliers at this plot point
	  for(unsigned i=0; i<Dim; i++)
	  {
	    outfile << lambda[i] << " ";
	  }

	  // output total pressure
	  outfile << p_total << " ";
	  
	  // output the traction at this point
	  for(unsigned i=0; i<Dim; i++)	
	  {
	    for(unsigned j=0; j<Dim; j++)	
	    {
	      outfile << stress_fe(i,j) << " ";
	    }
	  }

	  for(unsigned ising=0; ising<Nsingular_fct; ising++)
	  {
	    // Get gradient of scaled singular velocity functions
	    DenseMatrix<double> dudx_sing(Dim, Dim);
	    dudx_sing = Navier_stokes_sing_el_pt[ising]->gradient_of_singular_fct(x);

	    // compute the unscaled singular contribution to the strain-rate
	    DenseMatrix<double>strain_rate_sing(Dim, Dim, 0.0);
      
	    for (unsigned i=0; i<Dim; i++)
	    {
	      for(unsigned j=0; j<Dim; j++)
	      {
		strain_rate_sing(i,j) = 0.5*(dudx_sing(i,j) + dudx_sing(j,i));
	      }
	    }

	    // get contribution of singular pressure and singular velocity gradients to stress tensor
	    DenseMatrix<double> stress_sing(Dim, Dim);
	    stress_sing = (*bulk_el_pt->stress_fct_pt())(strain_rate_sing, p_sing[ising]);
	  
	    // output the traction at this point
	    for(unsigned i=0; i<Dim; i++)	
	    {
	      for(unsigned j=0; j<Dim; j++)	
	      {
		outfile << stress_sing(i,j) << " ";
	      }
	    }
	  }
	    outfile << std::endl;	  
	}
      }
      
    protected:

      /// \short Function to compute the shape and test functions and to return 
      /// the Jacobian of mapping between local and global (Eulerian)
      /// coordinates
      inline double shape_and_test(const Vector<double> &s, Shape &psi, Shape &test)
	const
      {
	//Find number of nodes
	unsigned n_node = nnode();

	//Get the shape functions
	shape(s,psi);

	//Set the test functions to be the same as the shape functions
	for(unsigned i=0; i<n_node; i++)
	{
	  test[i] = psi[i];
	}

	//Return the value of the jacobian
	return J_eulerian(s);
      }


      /// \short Function to compute the shape and test functions and to return 
      /// the Jacobian of mapping between local and global (Eulerian)
      /// coordinates
      inline double shape_and_test_at_knot(const unsigned &ipt,
					   Shape &psi, Shape &test)
	const
      {
	//Find number of nodes
	unsigned n_node = nnode();

	//Get the shape functions
	shape_at_knot(ipt,psi);

	//Set the test functions to be the same as the shape functions
	for(unsigned i=0; i<n_node; i++)
	{
	  test[i] = psi[i];
	}

	//Return the value of the jacobian
	return J_eulerian_at_knot(ipt);
      }

      const unsigned& boundary_id()
      {
	return Boundary_id;
      }
      
    private:

      /// \short Add the element's contribution to its residual vector.
      /// flag=1(or 0): do (or don't) compute the contribution to the
      /// Jacobian as well. 
      void fill_in_generic_residual_contribution_navier_stokes_sing(
	Vector<double> &residuals, DenseMatrix<double> &jacobian, 
	const unsigned& flag);
  
      /// The spatial dimension of the problem
      unsigned Dim;

      /// ID of the boundary this face element sits on
      unsigned Boundary_id;
      
      /// The index at which the Stokes unknown is stored at the nodes
      unsigned P_index_nst;

      /// \short The index at which the Lagrange multiplier that enforces
      /// the Dirichlet BC is stored at the nodes
      Vector<unsigned> Lambda_index;

      /// Desired boundary values at nodes
      DenseMatrix<double> Nodal_boundary_value;

      /// The number of singular functions which compose the solution
      unsigned Nsingular_fct;
      
      /// \short Indices of external Data that store the values of the amplitudes of
      /// the singular functions
      Vector<unsigned> C_external_data_index;
  
      /// \short Indices of values (within external Data) that store the
      /// values of the amplitudes of the singular functions
      Vector<unsigned> C_external_data_value_index;
  
      /// \short Vector of pointers to elements that store pointers to singular fcts
      /// (and their gradients etc.) as well as amplitudes
      Vector<ScalableSingularityForNavierStokesElement<ELEMENT>*> Navier_stokes_sing_el_pt;
    }; 

  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  
  //===========================================================================
  /// Constructor, takes the pointer to the "bulk" element, the 
  /// index of the fixed local coordinate and its value represented
  /// by an integer indicating which face we're on.
  /// Optional final arg is the identifier for the new values created
  /// by this face element
  //===========================================================================
  template<class ELEMENT>
    NavierStokesWithSingularityBCFaceElement<ELEMENT>::
    NavierStokesWithSingularityBCFaceElement(FiniteElement* const& bulk_el_pt, 
					     const int& face_index, 
					     const unsigned& id) : 
  FaceGeometry<ELEMENT>(), FaceElement(), Nsingular_fct(0)
  { 

    // set the identifier for the boundary we're on
    Boundary_id = id;
    
    /* // Initialise singular element pointer */
    /* Navier_stokes_sing_el_pt = 0; */

    // Let the bulk element build the FaceElement, i.e. setup the pointers 
    // to its nodes (by referring to the appropriate nodes in the bulk
    // element), etc.
    bulk_el_pt->build_face_element(face_index, this);
 
#ifdef PARANOID
    {
      //Check that the element is not a refineable 3d element
      ELEMENT* elem_pt = dynamic_cast<ELEMENT*>(bulk_el_pt);
      //If it's three-d
      if(elem_pt->dim() == 3)
      {
	//Is it refineable
	RefineableElement* ref_el_pt = dynamic_cast<RefineableElement*>(elem_pt);
	if(ref_el_pt != 0)
	{
	  if (this->has_hanging_nodes())
	  {
	    throw OomphLibError(
	      "This face element will not work correctly if nodes are hanging\n",
	      OOMPH_CURRENT_FUNCTION,
	      OOMPH_EXCEPTION_LOCATION);
	  }
	}
      }
    }
#endif   

    // Extract the dimension of the problem from the dimension of 
    // the first node
    Dim = this->node_pt(0)->ndim();

    // Set up P_index_nst. Initialise to Dim, (since we have Dim velocity components indexed
    // from zero, followed by the pressure) which probably won't change
    // in most cases, oh well, the price we pay for generality
    P_index_nst = Dim;

    // Cast to the appropriate NavierStokesEquation so that we can
    // find the index at which the variable is stored
    // We assume that the dimension of the full problem is the same
    // as the dimension of the node, if this is not the case you will have
    // to write custom elements, sorry
    switch(Dim)
    {
      //One dimensional problem
      case 1:
      {
	NavierStokesEquations<1>* eqn_pt = 
	  dynamic_cast<NavierStokesEquations<1>*>(bulk_el_pt);
	
	//If the cast has failed die
	if(eqn_pt == 0)
	{
	  std::string error_string =
	    "Bulk element must inherit from NavierStokesEquations.";
	  error_string += 
	    "Nodes are one dimensional, but cannot cast the bulk element to\n";
	  error_string += "NavierStokesEquations<1>\n.";
	  error_string += 
	    "If you desire this functionality, you must implement it yourself\n";
       
	  throw OomphLibError(error_string,
			      OOMPH_CURRENT_FUNCTION,
			      OOMPH_EXCEPTION_LOCATION);
	}
	//Otherwise read out the value
	else
	{
	  //Read the index from the (cast) bulk element
	  P_index_nst = eqn_pt->p_nodal_index_nst();
	}
      }
      break;
    
      //Two dimensional problem
      case 2:
      {
	NavierStokesEquations<2>* eqn_pt = 
	  dynamic_cast<NavierStokesEquations<2>*>(bulk_el_pt);
	//If the cast has failed die
	if(eqn_pt == 0)
	{
	  std::string error_string =
	    "Bulk element must inherit from NavierStokesEquations.";
	  error_string += 
	    "Nodes are two dimensional, but cannot cast the bulk element to\n";
	  error_string += "NavierStokesEquations<2>\n.";
	  error_string += 
	    "If you desire this functionality, you must implement it yourself\n";
       
	  throw OomphLibError(error_string,
			      OOMPH_CURRENT_FUNCTION,
			      OOMPH_EXCEPTION_LOCATION);
	}
	else
	{
	  //Read the index from the (cast) bulk element.
	  P_index_nst = eqn_pt->p_nodal_index_nst();
	}
      }
      break;
    
      //Three dimensional problem
      case 3:
      {
	NavierStokesEquations<3>* eqn_pt = 
	  dynamic_cast<NavierStokesEquations<3>*>(bulk_el_pt);
	
	//If the cast has failed die
	if(eqn_pt == 0)
	{
	  std::string error_string =
	    "Bulk element must inherit from NavierStokesEquations.";
	  error_string += 
	    "Nodes are three dimensional, but cannot cast the bulk element to\n";
	  error_string += "NavierStokesEquations<3>\n.";
	  error_string += 
	    "If you desire this functionality, you must implement it yourself\n";
       
	  throw OomphLibError(error_string,
			      OOMPH_CURRENT_FUNCTION,
			      OOMPH_EXCEPTION_LOCATION);       
	}
	else
	{
	  //Read the index from the (cast) bulk element.
	  P_index_nst = eqn_pt->p_nodal_index_nst();
	}
      }
      break;

      //Any other case is an error
      default:
	std::ostringstream error_stream; 
	error_stream <<  "Dimension of node is " << Dim 
		     << ". It should be 1,2, or 3!" << std::endl;
     
	throw OomphLibError(error_stream.str(),
			    OOMPH_CURRENT_FUNCTION,
			    OOMPH_EXCEPTION_LOCATION);
	break;
    }

    // Where is the extra dof representing the Lagrange multiplier stored?
    // Initially store number of values stored right now
    unsigned nnod = nnode();

    // Make space for Dim Lagrange multipliers
    Vector<unsigned> n_additional_values(nnod, Dim);
    this->add_additional_values(n_additional_values, id);

  } // end NavierStokesWithSingularityBCFaceElement constructor

  //===========================================================================
  /// \short Compute the element's residual vector and the Jacobian matrix.
  /// Adds this boundary face element's contribution to the equation which
  /// determines the Lagrange multipliers, and adds the Lagrange multiplier
  /// contribution to the bulk equations
  //===========================================================================
  template<class ELEMENT>
    void NavierStokesWithSingularityBCFaceElement<ELEMENT>::
    fill_in_generic_residual_contribution_navier_stokes_sing(
      Vector<double> &residuals, DenseMatrix<double> &jacobian, 
      const unsigned& flag) 
  {

    //Find out how many nodes there are
    const unsigned n_node = nnode();
     
    //Set up memory for the shape and test functions
    Shape psi(n_node), test(n_node);
     
    //Set the value of Nintpt
    const unsigned n_intpt = integral_pt()->nweight();
     
    //Set the Vector to hold local coordinates
    Vector<double> s(Dim-1);

    //Loop over the integration points
    //--------------------------------
    for(unsigned ipt=0; ipt<n_intpt; ipt++)
    {       
      //Assign values of s
      for(unsigned i=0; i<(Dim-1); i++) 
      {
	s[i] = integral_pt()->knot(ipt, i);
      }
       
      //Get the integral weight
      double w = integral_pt()->weight(ipt);
       
      //Find the shape and test functions and return the Jacobian
      //of the mapping
      double J = shape_and_test(s, psi, test);
       
      //Premultiply the weights and the Jacobian
      double W = w*J;
       
      //Calculate stuff at integration point
      Vector<double> interpolated_x(Dim, 0.0);
      Vector<double> u_fe(Dim, 0.0);
      Vector<double> u_bc(Dim, 0.0);	
      Vector<double> lambda(Dim, 0.0);
	
      for(unsigned l=0; l<n_node; l++)
      {
	// grab a pointer to the current node
	Node* node_pt = this->node_pt(l);

	// get the map which gives the starting nodal index for
	// the Lagrange multipliers associated with each boundary ID
	std::map<unsigned, unsigned> first_index = *(
	  dynamic_cast<BoundaryNodeBase*>(node_pt)->
	  index_of_first_value_assigned_by_face_element_pt() );

	for(unsigned i=0; i<Dim; i++)
	{
	  // get the nodal index, accounting for the dimension offset
	  unsigned lambda_index = first_index[Boundary_id] + i;
	    
	  // get the nodal values of the FE solution and the boundary conditions
	  u_fe[i] += this->nodal_value(l,i)    * psi[l];
	  u_bc[i] += Nodal_boundary_value(l,i) * psi[l];

	  // get the interpolated Lagrange multipliers
	  lambda[i] += this->nodal_value(l, lambda_index) * psi[l];

	  // get the interpolated position
	  interpolated_x[i] += this->nodal_position(l,i) * psi[l];	    
	}
	  
      }

      // Stuff related to singular fct
      Vector<double> u_sing(Dim, 0.0);
      
      Vector<Vector<double> > u_sing_unscaled(Nsingular_fct);

      // loop over the singular functions
      for(unsigned ising=0; ising<Nsingular_fct; ising++)
      {	
	if (Navier_stokes_sing_el_pt[ising] != 0)
	{
	  // get the scaled and unscaled versions of this singular function
	  Vector<double> u_sing_i, u_sing_unscaled_i;
          u_sing_i          = Navier_stokes_sing_el_pt[ising]->singular_fct(interpolated_x);

	  u_sing_unscaled[ising].resize(Dim);
	  u_sing_unscaled[ising] = Navier_stokes_sing_el_pt[ising]->unscaled_singular_fct(interpolated_x);

	  // add them onto the total
	  for(unsigned i=0; i<Dim; i++)
	  {
	    u_sing[i] += u_sing_i[i];	  
	  }
	}
      }
      //Now add to the appropriate equations

      // number of local equation which determines the singular amplitude
      int local_eqn_c = -1;

      // QUEHACERES think about how to loop this properly when we come to implementing the
      // analytic jacobian
      /* if (Navier_stokes_sing_el_pt != 0) */
      /* { */
      /* 	local_eqn_c = external_local_eqn(C_external_data_index, */
      /* 					 C_external_data_value_index); */
      /* } */

      //Loop over the test functions
      for(unsigned l=0; l<n_node; l++)
      {
	// grab a pointer to the current node
	Node* node_pt = this->node_pt(l);
	  
	// get the map which gives the starting nodal index for
	// the Lagrange multipliers associated with each boundary ID
	std::map<unsigned, unsigned> first_index = *(
	  dynamic_cast<BoundaryNodeBase*>(node_pt)->
	  index_of_first_value_assigned_by_face_element_pt() );

	// loop over the coordinate directions
	for(unsigned d=0; d<Dim; d++)
	{
	  // get the nodal index of the Lagrange multiplier for this
	  // coordinate direction and boundary ID
	  unsigned lambda_index = first_index[Boundary_id] + d;

	  // get the local Lagrange multiplier equation number 
	  int local_eqn_lagr = nodal_local_eqn(l, lambda_index);
	      
	  // QUEHACERES get this nodal index systematically, don't assume it starts at 0
	  int local_eqn_u_fe = nodal_local_eqn(l, d);

#ifdef PARANOID
	  // Lagrange multiplier active but u_fe pinned won't work!
	  if ( (local_eqn_lagr >= 0) && (local_eqn_u_fe < 0) )
	  {
	    throw OomphLibError(
	      "Lagrange multiplier active but u_fe pinned won't work!",
	      OOMPH_CURRENT_FUNCTION,
	      OOMPH_EXCEPTION_LOCATION);
	  }
#endif	    
	    	      
	  // Lagrange multiplier for BC residual: It's determined by enforcing
	  // that u_fe + C u_sing = u_bc
	  if(local_eqn_lagr >= 0)
	  {
	    // QUEHACERES don't fuck about with this sign, it's right!
	    residuals[local_eqn_lagr] -= (u_fe[d] + u_sing[d] - u_bc[d]) * test[l]*W;
	      
	    // Jacobian?
	    if (flag == 1)
	    {
	      for(unsigned l2=0; l2<n_node; l2++)
	      {
		// QUEHACERES again, get this index more systematically
		int local_unknown_u_fe = nodal_local_eqn(l2, d);
		if (local_unknown_u_fe >= 0)
		{
		  jacobian(local_eqn_lagr, local_unknown_u_fe) += psi[l2] * test[l]*W;
		}
	      }

	      // QUEHACERES needs to be looped properly when analytic jacobian is implemented
	      /* // Deriv. w.r.t. amplitude is simply the unscaled singular fct. */
	      /* if (local_eqn_c >= 0) */
	      /* { */
	      /* 	jacobian(local_eqn_lagr, local_eqn_c) += u_sing_unscaled[d] * test[l]*W; */
	      /* } */
	    }
	  }
         
	  // Contribution of Lagrange multiplier to bulk eqn:
	  if (local_eqn_u_fe >= 0)
	  {
	    residuals[local_eqn_u_fe] -= lambda[d] * test[l] * W;

	    // QUEHACERES need to review this code, never been tested
	    if (flag == 1)
	    {
	      for(unsigned l2=0; l2<n_node; l2++)
	      {
		// grab a pointer to the second node
		Node* node2_pt = this->node_pt(l2);

		// get the map which gives the starting nodal index for
		// the Lagrange multipliers associated with each boundary ID
		std::map<unsigned, unsigned> first_index2 = *(
		  dynamic_cast<BoundaryNodeBase*>(node2_pt)->
		  index_of_first_value_assigned_by_face_element_pt() );
		  
		// get the index of the Lagrange multiplier of the second node
		// associated with this face ID and direction 
		unsigned lambda_index2 = first_index2[Boundary_id] + d;
		int local_unknown_lambda = nodal_local_eqn(l2, lambda_index2);
		      
		if (local_unknown_lambda >= 0)
		{
		  jacobian(local_eqn_u_fe, local_unknown_lambda) += psi[l2] * test[l] * W;
		}
		    
	      }
	    }
	  }	    
	} // end loop over directions

      } // end loop over nodes
    } // end loop over integration points
  }

  
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////

  // hierher really need to tidy this up! Should only need one class 
  // for T and Q
  //
  //====================================================================
  /// New class. Mainly overloads output-related functions to add
  /// "singular function" (which is assumed to satisfy the Laplace
  /// equation; therefore no change to the governing (bulk) equations) 
  /// to the FE solution. 
  //====================================================================
  template<unsigned DIM, unsigned NNODE_1D>
    class MyTNavierStokesElement : public virtual TTaylorHoodElement<DIM>
  {
    
  public:

    typedef void (*ExactNonSingularFctPt)
      (const Vector<double>& x, Vector<double>& u, DenseMatrix<double>& grad_u);

    // function pointer for helper function which computes the stress
    typedef DenseMatrix<double> (*StressFctPt)(const DenseMatrix<double>& du_dx, const double& p);
        
    ExactNonSingularFctPt& exact_non_singular_fct_pt()
    {
      return Exact_non_singular_fct_pt;
    }  

    StressFctPt& stress_fct_pt()
    {
      return Stress_fct_pt;
    }
    
    /// Constructor
  MyTNavierStokesElement() : Nsingular_fct(0), Exact_non_singular_fct_pt(0)
    {
    }

    /// \short Return FE representation of function value u_navier_stokes(s) 
    /// plus scaled singular fct (if provided) at local coordinate s
    inline Vector<double> interpolated_u_total_navier_stokes(const Vector<double> &s) const
    {
      // FE part of the solution
      Vector<double> u_fe(DIM);
      
      // get the interpolated FE velocities
      TTaylorHoodElement<DIM>::interpolated_u_nst(s, u_fe);

      // get the interpolated FE pressure
      double p = TTaylorHoodElement<DIM>::interpolated_p_nst(s);

      // add pressure to the solution vector
      u_fe.push_back(p);

      for(unsigned ising=0; ising<Nsingular_fct; ising++)
      {
	// check if we're subtracting the singularity or not
	if (Navier_stokes_sing_el_pt[ising] != 0)
	{
	  // singular part of the solution
	  Vector<double> u_sing(DIM);

	  // interpolate the position
	  Vector<double> x(DIM);
	
	  for(unsigned i=0; i<DIM; i++)  
	  { 
	    x[i] = this->interpolated_x(s,i); 
	  }
	
	  // get singular part
	  u_sing = Navier_stokes_sing_el_pt[ising]->singular_fct(x);

	  // add singular part of the solution to the FE part to give the total
	  // computed solution
	  for(unsigned i=0; i<DIM+1; i++)
	  {
	    u_fe[i] += u_sing[i];
	  }
	}
      }
      
      return u_fe;
    } 

    /// \short Return FE representation of function value u_fe
    inline Vector<double> interpolated_u_fe_navier_stokes(const Vector<double>& s) const
    {
      // FE solution vector
      Vector<double> u_fe(DIM);

      // get FE velocity
      TTaylorHoodElement<DIM>::interpolated_u_nst(s, u_fe);

      // get FE pressure
      double p_fe = TTaylorHoodElement<DIM>::interpolated_p_nst(s);

      // add pressure to the solution vector
      u_fe.push_back(p_fe);
      
      return u_fe;
    } 

    void output_with_various_contributions(std::ostream& outfile, 
					   const Vector<double>& s)
    {
      Vector<double> x(DIM);
      for(unsigned i=0; i<DIM; i++) 
      {
	x[i] = this->interpolated_x(s,i);	
      }

      // regular part of the solution
      Vector<double> u_exact_non_sing(DIM+1, 0.0);

      DenseMatrix<double> dudx(DIM, DIM);

      // Overwrite with exact version!
      if (Exact_non_singular_fct_pt != 0)
      {
	Exact_non_singular_fct_pt(x, u_exact_non_sing, dudx);
      }
		
      // get the regular FE solution, and the full computed solution u = u_FE + u_sing
      Vector<double> u_fe(DIM+1, 0.0);
      Vector<double> u_fe_plus_sing(DIM+1, 0.0);

      u_fe           = this->interpolated_u_fe_navier_stokes(s);
      u_fe_plus_sing = this->interpolated_u_total_navier_stokes(s);


      // ==========================================
      // output total solution and FE bits
      // ==========================================
      // coordinates
      for(unsigned i=0; i<DIM; i++) 
      {
	outfile << x[i] << " ";
      }
      
      // output the total solution
      for(unsigned i=0; i<DIM+1; i++)
      {
	outfile << u_fe_plus_sing[i] << " ";
      }	
      // output the FE bit
      for(unsigned i=0; i<DIM+1; i++)
      {
	outfile << u_fe[i] << " ";
      }
      // ==========================================
            

      // ==========================================
      // output the singular bits
      // ==========================================

      for(unsigned ising=0; ising<Nsingular_fct; ising++)
      {
	// singular part of the solution
	Vector<double> u_sing(DIM+1, 0.0);

	// check we've got a singular function
	if (Navier_stokes_sing_el_pt[ising] != 0) 
	{ 
	  u_sing = Navier_stokes_sing_el_pt[ising]->singular_fct(x); 
	}
	
	for(unsigned i=0; i<DIM+1; i++)
	{
	  outfile << u_sing[i] << " ";
	}
      }
      // QUEHACERES leave stress out for the time being
      // FE stress and strain-rate tensors
      /* DenseMatrix<double> stress_fe(DIM, DIM); */
      /* DenseMatrix<double> strain_rate_fe(DIM, DIM); */

      /* // singular stress and strain-rate tensors */
      /* DenseMatrix<double> stress_sing(DIM, DIM); */
      /* DenseMatrix<double> stress_total(DIM, DIM); */
	
      /* DenseMatrix<double> strain_rate_sing(DIM, DIM); */
      /* DenseMatrix<double> strain_rate_total(DIM, DIM); */
	
      /* // get the strain rates */
      /* this->strain_rate(s, strain_rate_fe); */
	
      /* // Get gradient of scaled singular velocity functions       */
      /* DenseMatrix<double> dudx_sing(DIM, DIM); */
      /* if(Navier_stokes_sing_el_pt != 0) */
      /* { */
      /* 	dudx_sing = Navier_stokes_sing_el_pt->gradient_of_singular_fct(x); */
      /* } */
      
      /* // compute the unscaled singular contribution to the strain-rate       */
      /* for (unsigned i=0; i<DIM; i++) */
      /* { */
      /* 	for(unsigned j=0; j<DIM; j++) */
      /* 	{ */
      /* 	  strain_rate_sing(i,j)  = 0.5*(dudx_sing(i,j) + dudx_sing(j,i)); */
      /* 	  strain_rate_total(i,j) = strain_rate_fe(i,j) + strain_rate_sing(i,j); */
      /* 	} */
      /* } */
		
      /* // extract pressure from solution */
      /* double p_sing  = u_sing[DIM]; */
      /* double p_fe    = u_fe[DIM];	 */
      /* double p_total = u_fe_plus_sing[DIM]; */
	
      /* // compute stress from constitutive equation */
      /* stress_sing  = (this->stress_fct_pt())(strain_rate_sing, p_sing); */
      /* stress_fe    = (this->stress_fct_pt())(strain_rate_fe, p_fe); */
      /* stress_total = (this->stress_fct_pt())(strain_rate_total, p_total); */
      
      /* for(unsigned i=0; i<DIM; i++) */
      /* { */
      /* 	for(unsigned j=0; j<DIM; j++) */
      /* 	{ */
      /* 	  outfile << stress_total(i,j) << " "; */
      /* 	}  */
      /* } */
      /* for(unsigned i=0; i<DIM; i++) */
      /* { */
      /* 	for(unsigned j=0; j<DIM; j++) */
      /* 	{ */
      /* 	  outfile << stress_fe(i,j) << " "; */
      /* 	}  */
      /* } */
      /* for(unsigned i=0; i<DIM; i++) */
      /* { */
      /* 	for(unsigned j=0; j<DIM; j++) */
      /* 	{ */
      /* 	  outfile << stress_sing(i,j) << " "; */
      /* 	}  */
      /* } */
      outfile << std::endl;
    }
    
    /// Output with various contributions
    void output_with_various_contributions(std::ostream& outfile, 
					   const unsigned& nplot)
    {
      //Vector of local coordinates
      Vector<double> s(DIM);
   
      // Tecplot header info
      outfile << this->tecplot_zone_string(nplot);
   
      // Loop over plot points
      unsigned num_plot_points = this->nplot_points(nplot);
      for (unsigned iplot=0; iplot < num_plot_points; iplot++)
      {
	// Get local coordinates of plot point
	this->get_s_plot(iplot, nplot, s);

	// do the output
	output_with_various_contributions(outfile, s);	
      }
      outfile << std::endl;
      
      // Write tecplot footer (e.g. FE connectivity lists)
      this->write_tecplot_zone_footer(outfile, nplot);   
    }
    
    /// Pointer to element that stores ith singular fct 
    TemplateFreeScalableSingularityForNavierStokesElement*& navier_stokes_sing_el_pt(unsigned i) 
    { 
      return Navier_stokes_sing_el_pt[i]; 
    } 

    void add_singular_fct_pt(TemplateFreeScalableSingularityForNavierStokesElement* new_sing_el_pt)
    {
      // add this pointer to the list
      Navier_stokes_sing_el_pt.push_back(new_sing_el_pt);

      // increment the counter of the number of singular functions
      Nsingular_fct++;
    }

    unsigned nsingular_fct() const
    {
      return Nsingular_fct;
    }
      
  private:

    /// Number of singular functions which compose the total solution
    unsigned Nsingular_fct;
    
    /// Pointers to elements that stores singular fcts
    Vector<TemplateFreeScalableSingularityForNavierStokesElement*> Navier_stokes_sing_el_pt;

    /* TemplateFreeScalableSingularityForNavierStokesElement* Navier_stokes_sing_el_pt; */
    
    /// Pointer to exact non-singular fct (only for post-processing!)
    ExactNonSingularFctPt Exact_non_singular_fct_pt;

    /// Pointer to function which computes the stress
    StressFctPt Stress_fct_pt;
  };


  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////


  //=======================================================================
  /// Face geometry for the MyTNavierStokesElement elements: The spatial 
  /// dimension of the face elements is one lower than that of the
  /// bulk element but they have the same number of points
  /// along their 1D edges.
  //=======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
    class FaceGeometry<MyTNavierStokesElement<DIM,NNODE_1D> >: 
    public virtual TElement<DIM-1,NNODE_1D>
  {

  public:
 
    /// \short Constructor: Call the constructor for the
    /// appropriate lower-dimensional TElement
  FaceGeometry() : TElement<DIM-1,NNODE_1D>() {}

  };

  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////


  //=======================================================================
  /// Face geometry for the 1D MyTNavierStokesElement elements: Point elements
  //=======================================================================
  template<unsigned NNODE_1D>
    class FaceGeometry<MyTNavierStokesElement<1,NNODE_1D> >: 
    public virtual PointElement
    {

    public:
 
      /// \short Constructor: Call the constructor for the
      /// appropriate lower-dimensional TElement
    FaceGeometry() : PointElement() {} 

    };
  
}

#endif
