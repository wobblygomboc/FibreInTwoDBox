#ifndef OOMPH_PLATE_AS_GEOM_OBJ_HEADER
#define OOMPH_PLATE_AS_GEOM_OBJ_HEADER

//=====start_of_plate======================================================
/// \short Generic class which defines a plate as a geometric object with 
/// an intersection point with a high-resolution region around its ends
//=========================================================================
class Plate : public GeomObject
{
public:
  
  /// Empty constructor
Plate() : GeomObject(1,2) {}
  
  /// \short Position Vector at Lagrangian coordinate zeta 
  virtual void position(const Vector<double>& zeta, Vector<double>& r) const = 0; 
    
  /// \short Position Vector at Lagrangian coordinate zeta  at time level t
  /// (t=0: present; t>0: previous level). 
  virtual void position(const unsigned& t, const Vector<double>& zeta,
			Vector<double>& r) const = 0;
  
  /// \short Position of the intersection point between the plate and the
  /// circle bounding the high-resolution region on the left side
  void hi_res_intersection_position_start(Vector<double>& x)
  {
    Vector<double> zeta(1);
    zeta[0] = Hi_res_zeta_start;
    position(zeta, x);
  }

  /// \short Position of the intersection point between the plate and the
  /// circle bounding the high-resolution region on the right side  
  void hi_res_intersection_position_end(Vector<double>& x)
  {
    Vector<double> zeta(1);
    zeta[0] = Hi_res_zeta_end;
    position(zeta, x);
  }

  double zeta_start()
  {
    return Zeta_start;
  }

  double zeta_end()
  {
    return Zeta_end;
  }

  double hi_res_zeta_start()
  {
    return Hi_res_zeta_start;
  }

  double hi_res_zeta_end()
  {
    return Hi_res_zeta_end;
  }

  unsigned nsegment()
  {
    return Nsegment;
  }
  
protected:
 

  /// Zeta for the start and end of the curve
  double Zeta_start;
  double Zeta_end;
  
  /// Intrinsic coordinate of intersection with high-resolution region
  double Hi_res_zeta_start;
  double Hi_res_zeta_end;

  // number of segments to use when meshing
  unsigned Nsegment;
  
}; // end of Plate class

//=====start_of_circular_plate============================================
/// \short Simple circle in 2D space.
/// \f[ x = X_c + R \cos(\zeta)  \f]
/// \f[ y = Y_c + R \sin(\zeta)  \f]
//=========================================================================
class CircularPlate : public Plate
{
public:
  
  /// Constructor:  Pass x and y-coords of centre and radius
CircularPlate(const double& x_c, const double& y_c, const double& r,
	      const double& zeta_start, const double& zeta_end,
	      double hi_res_zeta = 0) :
  X_c(x_c), Y_c(y_c), R(r)
  {
    // set the intrinsic coordinates for the start and end of the plate
    Zeta_start = zeta_start;
    Zeta_end   = zeta_end;
      
    // set the intrinsic coordinates of the intersection points of the hi-res regions
    Hi_res_zeta_start = Zeta_start + hi_res_zeta;
    Hi_res_zeta_end   = Zeta_end   - hi_res_zeta;

    Nsegment = 20;
  }
  
  /// \short Position Vector at Lagrangian coordinate zeta 
  void position(const Vector<double>& zeta, Vector<double>& r) const
  {
    // Position vector
    r[0] = X_c + R*cos(zeta[0]);
    r[1] = Y_c + R*sin(zeta[0]);
  }
  
  /// \short Position Vector at Lagrangian coordinate zeta  at time level t
  /// (t=0: present; t>0: previous level). Steady object, so we 
  /// simply forward the call to the steady version.
  void position(const unsigned& t, const Vector<double>& zeta,
		Vector<double>& r) const
  {
    position(zeta,r);
  }

protected:
  
  /// X-coordinate of centre
  double X_c;
  
  /// Y-coordinate of centre
  double Y_c;
  
  /// Radius
  double R;
  
}; // end of Plate class

class FlatPlate : public Plate
{
public:
  
FlatPlate(const Vector<double>& start_coords, const Vector<double>& end_coords,
	  double hi_res_zeta = 0) : Start_coords(start_coords)
  {   
    // get the length of the plate
    double length = sqrt(pow(end_coords[0] - start_coords[0], 2) +
			 pow(end_coords[1] - start_coords[1], 2));

    Zeta_start = 0;
    Zeta_end   = length;

    // set the intrinsic coordinates of the intersection points of the hi-res regions
    Hi_res_zeta_start = Zeta_start + hi_res_zeta;
    Hi_res_zeta_end   = Zeta_end   - hi_res_zeta;
    
    Unit_vector.resize(2);
    Unit_vector[0] = (end_coords[0] - start_coords[0]) / length;
    Unit_vector[1] = (end_coords[1] - start_coords[1]) / length;

    Nsegment = 5;
  }

  void position(const Vector<double>& zeta, Vector<double>& r) const
  {
    // position is simply how far along we are
    r[0] = Start_coords[0] + zeta[0] * Unit_vector[0];
    r[1] = Start_coords[1] + zeta[0] * Unit_vector[1];
  }  

  /// \short Position Vector at Lagrangian coordinate zeta  at time level t
  /// (t=0: present; t>0: previous level). Steady object, so we 
  /// simply forward the call to the steady version.
  void position(const unsigned& t, const Vector<double>& zeta,
		Vector<double>& r) const
  {
    position(zeta,r);
  }
  

protected:

  // Eulerian coordinates of the start and end of the plate
  Vector<double> Start_coords;

  // unit vector tangent to the plate  
  Vector<double> Unit_vector;
};


#endif
