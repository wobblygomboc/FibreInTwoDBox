#ifndef OOMPH_UNSTRUCTURED_MOFFAT_MESH_HEADER
#define OOMPH_UNSTRUCTURED_MOFFAT_MESH_HEADER
// Christian Vaquero-Stainer
// 21/11/19
// 
// domain is a rectangular box with x \in [L_edge_x, R_edge_x] and
// y \in [bottom_edge_y, top_edge_y]
// with a curved plate (B3) inside running from (-plate_radius,0) to (+plate_radius,0).
//
// Two regions are defined, R0 is y<=0 and R1 is y<0
// which allows for a pressure jump across the internal boundary.
//
// ==================================================================
// boundaries (B) and regions (R) are labelled as follows:
//
//     (L_edge_x, top_edge_y)       (R_edge_x, top_edge_y) 
//                     V_______________V
//                     |               |
//                  B0 |      R1       | 
//                     |               |
//       (L_edge_x, 0) |.....     .....| (R_edge_x,0) 
//                     , B2  \___/ B4  , 
//                     |      B3       |
//                     |               |
//                     ,      R0       ,
//                     |_,_,_,_,_,_,_,_|
//  (L_edge_x, bottom_edge_y)     B1    (R_edge_x, bottom_edge_y)

       
// sign function
template <typename T>
int sgn(T val)
{
  return (T(0) < val) - (val < T(0));
}

// return unsigned double
double udouble(double val)
{
  return val * sgn(val);
}

//=====================================================================
/// Helper function to build the mesh; assumed to live in namespace
/// where dimensions of mesh are defined
//=====================================================================
/* TriangleMesh<ELEMENT>* build_the_mesh(const double& uniform_element_area) */
template<class ELEMENT>
RefineableTriangleMesh<ELEMENT>* build_the_mesh(const double& uniform_element_area,
						Plate* plate_pt)
{
 
  // The domain is defined by 2 outer boundaries
  // (the upper and lower sections), and two internal boundaries,
  // one representing the plate and one to complete the region boundary
  Vector<TriangleMeshCurveSection*> outer_boundary_polyline_pt(2);
  
  // =========================================================
  // Boundary 0: Upper outer boundary

  // define start and end points of straight line section
  Vector <Vector <double> > bound_0(4, Vector<double>(2));
  
  bound_0[0][0] = L_edge_x;
  bound_0[0][1] = 0.0;
  
  bound_0[1][0] = L_edge_x;
  bound_0[1][1] = top_edge_y;
  
  bound_0[2][0] = R_edge_x;
  bound_0[2][1] = top_edge_y;
  
  bound_0[3][0] = R_edge_x;
  bound_0[3][1] = 0.0;

  // create the polygonal line for the upper half of the outer boundary
  TriangleMeshPolyLine* upper_boundary_pt =
    new TriangleMeshPolyLine(bound_0, Outer_boundary_upper);

  // add it to the array of outer boundaries
  outer_boundary_polyline_pt[0] = upper_boundary_pt;

  // =========================================================
  // Boundary 1: Lower outer boundary

  // define start and end points of straight line section
  Vector <Vector <double> > bound_1(4, Vector<double>(2));
  
  bound_1[0][0] = R_edge_x;
  bound_1[0][1] = 0.0;
  
  bound_1[1][0] = R_edge_x;
  bound_1[1][1] = bottom_edge_y;
  
  bound_1[2][0] = L_edge_x;
  bound_1[2][1] = bottom_edge_y;
  
  bound_1[3][0] = L_edge_x;
  bound_1[3][1] = 0.0;
  
  outer_boundary_polyline_pt[1] = new TriangleMeshPolyLine(bound_1, Outer_boundary_lower);  

  // =========================================================
  // Inner boundaries
  unsigned n_inner_boundaries = include_hi_res_regions ? 7 : 3;
  Vector<TriangleMeshOpenCurve*> inner_open_boundary_pt(n_inner_boundaries);

  // get the coordinates of the intersection between the plate and the left high-res bounding
  // circle
  Vector<double> r_intersection_left(2);
  Vector<double> r_intersection_right(2);
  plate_pt->hi_res_intersection_position_start(r_intersection_left);
  plate_pt->hi_res_intersection_position_end(r_intersection_right);
  
  // compute the radius of the high-resolution region from the intrinsic coordinate
  // of the intersection 
  double high_res_region_radius_left  =
    sqrt(pow(-plate_radius - r_intersection_left[0],2) + pow(r_intersection_left[1],2));
  double high_res_region_radius_right =
    sqrt(pow(plate_radius - r_intersection_right[0],2) + pow(r_intersection_right[1],2));
  
  //  Boundary 2: left region boundary
  
  // define start and end points of straight line section
  Vector < Vector <double> > bound_2(3, Vector<double>(2));
 
  bound_2[0][0] = L_edge_x;
  bound_2[0][1] = 0.0;

  // extra node, to allow connections with the circle defining the high-res
  // region around the end of the plate
  bound_2[1][0] = -plate_radius - high_res_region_radius_left;
  bound_2[1][1] = 0.0;
  
  bound_2[2][0] = -plate_radius;
  bound_2[2][1] = 0.0;

  TriangleMeshPolyLine* left_region_boundary_polyline_pt =
    new TriangleMeshPolyLine(bound_2, Inner_region_boundary_left);
  
  // Do the connection with the outer (upper) boundary
  left_region_boundary_polyline_pt->connect_initial_vertex_to_polyline(upper_boundary_pt, 0);

  Vector<TriangleMeshCurveSection*> inner_boundary_line1_pt(1);
  inner_boundary_line1_pt[0] = left_region_boundary_polyline_pt;
  
  // add to array of inner boundaries
  inner_open_boundary_pt[0] = new TriangleMeshOpenCurve(inner_boundary_line1_pt);

  // =========================================================
  // Boundary 3: The plate
    
  // Number of segments used for representing the curvlinear internal boundary
  unsigned n_segments = plate_pt->nsegment();
 
  TriangleMeshCurviLine* plate_boundary_curviline_pt =
    new TriangleMeshCurviLine(plate_pt, plate_pt->zeta_start(),
			      plate_pt->zeta_end(), n_segments, Inner_plate_boundary);

  // QUEHACERES
  double tolerance = 1e-14; // 2e-13;
  
  // Do the connection with the right vertex of the left region boundary
  plate_boundary_curviline_pt->
    connect_initial_vertex_to_polyline(left_region_boundary_polyline_pt, 2, tolerance);

  Vector<TriangleMeshCurveSection*> inner_boundary_line2_pt(1);
  inner_boundary_line2_pt[0] = plate_boundary_curviline_pt;
  
  // add to array of inner boundaries
  inner_open_boundary_pt[1] = new TriangleMeshOpenCurve(inner_boundary_line2_pt);
  
  // =========================================================
  // Hi-res left region  
  Vector<TriangleMeshCurveSection*> hi_res_left_boundary_curviline_pt(2);
      
  // =========================================================
  // Boundary 4: Right region boundary
 
  // define start and end points of straight line section
  Vector < Vector <double> > bound_4(3, Vector<double>(2));
    
  bound_4[0][0] = plate_radius;
  bound_4[0][1] = 0.0;

  // extra node, to allow connections with the circle defining the high-res
  // region around the end of the plate
  bound_4[1][0] = plate_radius + high_res_region_radius_right;
  bound_4[1][1] = 0.0;

  bound_4[2][0] = R_edge_x;
  bound_4[2][1] = 0.0;
  
  TriangleMeshPolyLine* right_region_boundary_polyline_pt =
    new TriangleMeshPolyLine(bound_4, Inner_region_boundary_right);
  
  // connect start to the plate and the end to the outer (upper) boundary
  right_region_boundary_polyline_pt ->
    connect_initial_vertex_to_curviline(plate_boundary_curviline_pt,
					plate_pt->zeta_end(), tolerance);

  // State the vertex number for connection on the destination boundaries
  unsigned vertex_to_connect_final = 3;
  
  right_region_boundary_polyline_pt ->
    connect_final_vertex_to_polyline(upper_boundary_pt, vertex_to_connect_final);

  Vector<TriangleMeshCurveSection*> inner_boundary_line3_pt(1);
  inner_boundary_line3_pt[0] = right_region_boundary_polyline_pt;
  
  // add to array of inner boundaries
  inner_open_boundary_pt[2] = new TriangleMeshOpenCurve(inner_boundary_line3_pt);

  // =========================================================
  if(include_hi_res_regions)
  {
    // Boundary 5: Upper half of left hi-res region

    // centre of the circle
    double left_edge_x = -plate_radius;
    double left_edge_y = 0.0;

    Circle* left_high_res_circle_pt =
      new Circle(left_edge_x, left_edge_y, high_res_region_radius_left);

    // The intrinsic coordinates for the beginning and end of the curve
    double s_hi_res_left_start = atan2(r_intersection_left[1] - left_edge_y,
				       r_intersection_left[0] - left_edge_x);
  
    double s_hi_res_left_end   = MathematicalConstants::Pi;
  
    // build the upper semi-circle
    n_segments = 20;
    TriangleMeshCurviLine* high_res_left_upper_curviline_pt =
      new TriangleMeshCurviLine(left_high_res_circle_pt,
				s_hi_res_left_start, s_hi_res_left_end, n_segments,
				Inner_hi_res_region_left_upper);

    // add to array of closed curves
    hi_res_left_boundary_curviline_pt[0] = high_res_left_upper_curviline_pt;

    // connect the start and end points to the plate and the left internal region boundary
    high_res_left_upper_curviline_pt->
      connect_initial_vertex_to_curviline(plate_boundary_curviline_pt,
					  plate_pt->hi_res_zeta_start(), tolerance);
    high_res_left_upper_curviline_pt->
      connect_final_vertex_to_polyline(left_region_boundary_polyline_pt, 1);

    Vector<TriangleMeshCurveSection*> high_res_curve_section_pt(1);
    high_res_curve_section_pt[0] = high_res_left_upper_curviline_pt;
  
  
    inner_open_boundary_pt[3] = new TriangleMeshOpenCurve(high_res_curve_section_pt);
  
  
    // =========================================================
    // Boundary 6: Lower half of left hi-res region
  
    // build the lower semi-circle
    TriangleMeshCurviLine* high_res_left_lower_curviline_pt =
      new TriangleMeshCurviLine(left_high_res_circle_pt,
				-s_hi_res_left_end, s_hi_res_left_start, n_segments,
				Inner_hi_res_region_left_lower);

    // connect its end points
    high_res_left_lower_curviline_pt->connect_initial_vertex_to_curviline(
      high_res_left_upper_curviline_pt, s_hi_res_left_end);
    
    high_res_left_lower_curviline_pt->connect_final_vertex_to_curviline(
      high_res_left_upper_curviline_pt, s_hi_res_left_start);

  
    hi_res_left_boundary_curviline_pt[1] = high_res_left_lower_curviline_pt;

    Vector<TriangleMeshCurveSection*> high_res_lower_curve_section_pt(1);
    high_res_lower_curve_section_pt[0] = high_res_left_lower_curviline_pt;
  
    if(include_hi_res_regions)
    { 
      inner_open_boundary_pt[4] = new TriangleMeshOpenCurve(high_res_lower_curve_section_pt);
    }
  
    // =========================================================
    // Hi-res right region
  
    Vector<TriangleMeshCurveSection*> hi_res_right_boundary_curviline_pt(3);

    // =========================================================
    // Boundary 7: Upper half of right hi-res region

    double right_edge_x = plate_radius;
    double right_edge_y = 0.0;

    Circle* right_high_res_circle_pt =
      new Circle(right_edge_x, right_edge_y, high_res_region_radius_right);

    // The intrinsic coordinates for the beginning and end of the curve
    double s_hi_res_right_start   = 0.0;
    double s_hi_res_right_end     = atan2(r_intersection_right[1], r_intersection_right[0] - plate_radius);
    double s_hi_res_right_end_2pi =
      Additional_Maths_Functions::atan2pi(r_intersection_right[1], r_intersection_right[0] - plate_radius);
  
    // build the upper semi-circle
    TriangleMeshCurviLine* high_res_right_upper_curviline_pt =
      new TriangleMeshCurviLine(right_high_res_circle_pt,
				s_hi_res_right_start, s_hi_res_right_end_2pi, n_segments,
				Inner_hi_res_region_right_upper);

    // add to array of closed curves
    hi_res_right_boundary_curviline_pt[0] = high_res_right_upper_curviline_pt;

    // connect the start and end points to the right internal region boundary and the plate
    high_res_right_upper_curviline_pt->
      connect_final_vertex_to_curviline(plate_boundary_curviline_pt,
					plate_pt->hi_res_zeta_end(), tolerance);
    high_res_right_upper_curviline_pt->
      connect_initial_vertex_to_polyline(right_region_boundary_polyline_pt, 1);

    // add it to the inner regions array
    Vector<TriangleMeshCurveSection*> hi_res_right_boundary_curve_sec_pt(1);
    hi_res_right_boundary_curve_sec_pt[0] = high_res_right_upper_curviline_pt;

    if(include_hi_res_regions)
    {
      inner_open_boundary_pt[5] = new TriangleMeshOpenCurve(hi_res_right_boundary_curve_sec_pt);
    }
  
    // =========================================================
    // Boundary 8: Lower half of right hi-res region

    // catch the case where the plate is either practically or exactly flat, and the
    // hi-res intersection point is positive
    if(s_hi_res_right_end > 0 && s_hi_res_right_end < MathematicalConstants::Pi + tolerance)
    {
      s_hi_res_right_end = -s_hi_res_right_end;
    }
  
    // build the lower semi-circle
    TriangleMeshCurviLine* high_res_right_lower_curviline_pt =
      new TriangleMeshCurviLine(right_high_res_circle_pt,
				s_hi_res_right_end, s_hi_res_right_start, n_segments,
				Inner_hi_res_region_right_lower);
  
    // connect it's end points
    high_res_right_lower_curviline_pt->
      connect_final_vertex_to_curviline(high_res_right_upper_curviline_pt,
					s_hi_res_right_start, tolerance);

    high_res_right_lower_curviline_pt->
      connect_initial_vertex_to_curviline(high_res_right_upper_curviline_pt, // plate_boundary_curviline_pt,
					  s_hi_res_right_end_2pi, tolerance); // connection_angle

    hi_res_right_boundary_curviline_pt[2] = high_res_right_lower_curviline_pt;
  
    Vector<TriangleMeshCurveSection*> hi_res_right_lower_boundary_curve_sec_pt(1);
    hi_res_right_lower_boundary_curve_sec_pt[0] = high_res_right_lower_curviline_pt;

    inner_open_boundary_pt[6] = new TriangleMeshOpenCurve(hi_res_right_lower_boundary_curve_sec_pt);
    
  } // end if(include_hi_res_regions)
  
  // **********************************************************************************
  
  // =========================================================  
  // Point identifying regions above plate
  Vector<double> region1_point(2);

  // left hi-res region
  Vector<double> region2_point(2);
  Vector<double> region3_point(2);

  // right hi-res region
  Vector<double> region4_point(2);
  Vector<double> region5_point(2);
  
  // Define the coordinates of a point in the upper region
  region1_point[0] = 0.0;
  region1_point[1] = top_edge_y/2.0;

  // point in the upper half of the left high-res region
  region2_point[0] = -plate_radius;
  region2_point[1] = high_res_region_radius_left/2.0;

  // point in the lower half of the left high-res region
  region3_point[0] = -plate_radius - high_res_region_radius_left/2.0;
  region3_point[1] = -high_res_region_radius_left/2.0;

  // point in the upper half of the right high-res region
  region4_point[0] = plate_radius;
  region4_point[1] = high_res_region_radius_right/2.0;

  // point in the lower half of the right high-res region
  region5_point[0] = plate_radius + high_res_region_radius_right/2.0;
  region5_point[1] = -high_res_region_radius_left/2.0;

  // ----------------------------------------------------
  // Create the triangle mesh polygon for outer boundary
  // ----------------------------------------------------
  TriangleMeshPolygon *outer_polygon = new TriangleMeshPolygon(outer_boundary_polyline_pt);
 
  // Enable redistribution of polylines
  outer_polygon -> enable_redistribution_of_segments_between_polylines();
 
  // Set the pointer
  TriangleMeshClosedCurve* outer_boundary_closed_curve_pt = outer_polygon;

  // Now build the mesh
  //===================
 
  // Use the TriangleMeshParameters object for helping on the manage of the
  // TriangleMesh parameters
  TriangleMeshParameters triangle_mesh_parameters(outer_boundary_closed_curve_pt);
 
  // Specify the maximum area element
  triangle_mesh_parameters.element_area() = uniform_element_area;
 
  // Specify the internal open boundaries
  triangle_mesh_parameters.internal_open_curves_pt() = inner_open_boundary_pt;
 
  // Identify regions
  triangle_mesh_parameters.add_region_coordinates(Region_upper, region1_point);

  if(include_hi_res_regions)
  {
    triangle_mesh_parameters.add_region_coordinates(Region_hi_res_left_upper,  region2_point);
    triangle_mesh_parameters.add_region_coordinates(Region_hi_res_left_lower,  region3_point);
    triangle_mesh_parameters.add_region_coordinates(Region_hi_res_right_upper, region4_point);
    triangle_mesh_parameters.add_region_coordinates(Region_hi_res_right_lower, region5_point);

    // set the higher resolution for the regions around the plate ends
    triangle_mesh_parameters.set_target_area_for_region(Region_hi_res_left_upper,  hi_res_element_area);
    triangle_mesh_parameters.set_target_area_for_region(Region_hi_res_left_lower,  hi_res_element_area);
    triangle_mesh_parameters.set_target_area_for_region(Region_hi_res_right_upper, hi_res_element_area);
    triangle_mesh_parameters.set_target_area_for_region(Region_hi_res_right_lower, hi_res_element_area);
  }
  
  // Build the damn thing
  RefineableTriangleMesh<ELEMENT>* Bulk_mesh_pt =
    new RefineableTriangleMesh<ELEMENT>(triangle_mesh_parameters);
  /* TriangleMesh<ELEMENT>* Bulk_mesh_pt = */
  /*   new TriangleMesh<ELEMENT>(triangle_mesh_parameters); */
 
  return Bulk_mesh_pt;
}

#endif
