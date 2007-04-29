/*  
    Geodesy and Mapping C++ Library (GNU GaMa / GaMaLib)
    Copyright (C) 2007  Ales Cepek <cepek@gnu.org>

    This file is part of the GNU Gama C++ library
    
    This library is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this library; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*
 *  $Id: dataparser_g3adj.cpp,v 1.5 2007/04/29 15:02:34 cepek Exp $
 */



#include <gnu_gama/g3/g3_adjres.h>
#include <gnu_gama/xml/dataparser.h>
#include <gnu_gama/gon2deg.h>
#include <gnu_gama/radian.h>
#include <cstring>

using namespace std;
using namespace GNU_gama;

namespace GNU_gama {

  struct DataParser_g3adj {

    g3::AdjustmentResults* adj;

    DataParser_g3adj() : adj(0)
    {
    }
    ~DataParser_g3adj()
    {  
      delete adj;
    }

  };

}


void DataParser::close_g3adj()
{
  delete g3adj;
}

void DataParser::init_g3adj()
{
  g3adj = new DataParser_g3adj;


  // .....  <g3-adjustment-results>  .................................

  init(s_gama_data, t_g3_adj_results,
       s_g3a_adj_results, s_g3a_o_observations_end, 0,
       &DataParser::g3a_s_adj_results, 0, &DataParser::g3a_s_adj_results);

  // .....  <g3-adjustment-results>  <adjustment-statistics>  ........

  init(s_g3a_adj_results, t_adj_statistics,
       s_g3a_s_statistics, 0, s_g3a_s_statistics_end,
       0, 0, 0);

  init(s_g3a_s_statistics, t_algorithm,
       s_g3a_s_s_algorithm, 0, 0,
       0, &DataParser::add_text, &DataParser::g3a_s_algorithm);

  init(s_g3a_s_statistics, t_ellipsoid,
       s_g3a_s_ellipsoid, 0, 0,
       0, 0, 0);

  init(s_g3a_s_ellipsoid, t_caption,
       s_g3a_s_ellipsoid_cap, 0, 0,
       0, &DataParser::add_text, &DataParser::g3a_s_ell_caption);

  init(s_g3a_s_ellipsoid, t_id,
       s_g3a_s_ellipsoid_id, 0, 0,
       0, &DataParser::add_text, &DataParser::g3a_s_ell_id);

  init(s_g3a_s_ellipsoid, t_a,
       s_g3a_s_ellipsoid_a, 0, 0,
       0, &DataParser::add_text, &DataParser::g3a_s_ell_a);

  init(s_g3a_s_ellipsoid, t_b,
       s_g3a_s_ellipsoid_b, 0, 0,
       0, &DataParser::add_text, &DataParser::g3a_s_ell_b);

  init(s_g3a_s_statistics, t_parameters,
       s_g3a_s_parameters, 0, 0,
       0, &DataParser::add_text, &DataParser::g3a_s_parameters);

  init(s_g3a_s_statistics, t_equations,
       s_g3a_s_equations, 0, 0,
       0, &DataParser::add_text, &DataParser::g3a_s_equations);

  init(s_g3a_s_statistics, t_defect,
       s_g3a_s_defect, 0, 0,
       0, &DataParser::add_text, &DataParser::g3a_s_defect);

  init(s_g3a_s_statistics, t_redundancy,
       s_g3a_s_redundancy, 0, 0,
       0, &DataParser::add_text, &DataParser::g3a_s_redundancy);

  init(s_g3a_s_statistics, t_sum_of_squares,
       s_g3a_s_sum_of_squares, 0, 0,
       0, &DataParser::add_text, &DataParser::g3a_s_sum_of_squares);

  init(s_g3a_s_statistics, t_apriori_var,
       s_g3a_s_apriori_var, 0, 0,
       0, &DataParser::add_text, &DataParser::g3a_s_apriori_var);

  init(s_g3a_s_statistics, t_aposteriori_var,
       s_g3a_s_aposteriori_var, 0, 0,
       0, &DataParser::add_text, &DataParser::g3a_s_aposteriori_var);

  init(s_g3a_s_statistics, t_variance_factor,
       s_g3a_s_variance_factor, 0, 0,
       0, &DataParser::add_text, &DataParser::g3a_s_variance_factor);

  init(s_g3a_s_statistics, t_design_m_graph,
       s_g3a_s_design_m_graph, 0, 0,
       0, &DataParser::add_text, &DataParser::g3a_s_design_m_graph);

  // .....  <g3-adjustment-results>  <adjustment-results>  ........

  init(/*s_g3a_adj_results*/s_g3a_s_statistics_end, t_adj_results,
       s_g3a_r_adj_results, 0, s_g3a_r_adj_results_end,
       0, 0, 0);

  // .....  <g3-adjustment-results>  <adjustment-results>  <point>


  init(s_g3a_r_adj_results, t_point,
       s_g3a_r_point, s_g3a_r_point_after_u, 0,
       &DataParser::g3a_r_point, 0, &DataParser::g3a_r_point);

  init(s_g3a_r_point, t_id,
       s_g3a_r_point_id, 0, s_g3a_r_point_after_id,
       0,  &DataParser::add_text, &DataParser::g3a_r_point_id);

  // .....  <n> / <e> / <u> .......................................

  init(s_g3a_r_point_after_id, t_n,
       s_g3a_r_point_n, 0, s_g3a_r_point_after_n,
       0, 0, 0);

  init(s_g3a_r_point_after_n, t_e,
       s_g3a_r_point_e, 0, s_g3a_r_point_after_e,
       0, 0, 0);

  init(s_g3a_r_point_after_e, t_u,
       s_g3a_r_point_u, 0, s_g3a_r_point_after_u,
       0, 0, 0);

  // .....  <point> ... <n/e/u> <fixed/>  .........................

  init(s_g3a_r_point_n, t_fixed,
       s_g3a_r_point_n_type, 0, 0,
       0, 0, &DataParser::g3a_r_point_n_type);

  init(s_g3a_r_point_e, t_fixed,
       s_g3a_r_point_e_type, 0, 0,
       0, 0, &DataParser::g3a_r_point_e_type);

  init(s_g3a_r_point_u, t_fixed,
       s_g3a_r_point_u_type, 0, 0,
       0, 0, &DataParser::g3a_r_point_u_type);

  // .....  <point> ... <n/e/u> <free/>  ..........................

  init(s_g3a_r_point_n, t_free,
       s_g3a_r_point_n_type, 0, 0,
       0, 0, &DataParser::g3a_r_point_n_type);

  init(s_g3a_r_point_e, t_free,
       s_g3a_r_point_e_type, 0, 0,
       0, 0, &DataParser::g3a_r_point_e_type);

  init(s_g3a_r_point_u, t_free,
       s_g3a_r_point_u_type, 0, 0,
       0, 0, &DataParser::g3a_r_point_u_type);

  // .....  <point> ... <n/e/u> <constr/>  ........................

  init(s_g3a_r_point_n, t_constr,
       s_g3a_r_point_n_type, 0, 0,
       0, 0, &DataParser::g3a_r_point_n_type);

  init(s_g3a_r_point_e, t_constr,
       s_g3a_r_point_e_type, 0, 0,
       0, 0, &DataParser::g3a_r_point_e_type);

  init(s_g3a_r_point_u, t_constr,
       s_g3a_r_point_u_type, 0, 0,
       0, 0, &DataParser::g3a_r_point_u_type);

  // .....  <point> ... <n/e/u> <unused/>  ........................

  init(s_g3a_r_point_n, t_unused,
       s_g3a_r_point_n_type, 0, 0,
       0, 0, &DataParser::g3a_r_point_n_type);

  init(s_g3a_r_point_e, t_unused,
       s_g3a_r_point_e_type, 0, 0,
       0, 0, &DataParser::g3a_r_point_e_type);

  init(s_g3a_r_point_u, t_unused,
       s_g3a_r_point_u_type, 0, 0,
       0, 0, &DataParser::g3a_r_point_u_type);

  // ..... <point> ... <n/e/u> <dn/de/du/index>  ..................

  init(s_g3a_r_point_n, t_dn,
       s_g3a_r_point_n_dn, 0, 0,
       0, &DataParser::add_text, &DataParser::g3a_r_point_n_dn);

  init(s_g3a_r_point_e, t_de,
       s_g3a_r_point_e_de, 0, 0,
       0, &DataParser::add_text, &DataParser::g3a_r_point_e_de);

  init(s_g3a_r_point_u, t_du,
       s_g3a_r_point_u_du, 0, 0,
       0, &DataParser::add_text, &DataParser::g3a_r_point_u_du);

  init(s_g3a_r_point_n, t_ind,
       s_g3a_r_point_n_ind, 0, 0,
       0, &DataParser::add_text, &DataParser::g3a_r_point_n_ind);

  init(s_g3a_r_point_e, t_ind,
       s_g3a_r_point_e_ind, 0, 0,
       0, &DataParser::add_text, &DataParser::g3a_r_point_e_ind);

  init(s_g3a_r_point_u, t_ind,
       s_g3a_r_point_u_ind, 0, 0,
       0, &DataParser::add_text, &DataParser::g3a_r_point_u_ind);

  // ..... <point> <cnn>, <cne>, <cnu>, <cee>, <ceu>, <cuu> ........

  init(s_g3a_r_point_after_u, t_cnn,
       s_g3a_r_point_cnn, 0, 0,
       0, &DataParser::add_text, &DataParser::g3a_r_point_cnn);

  init(s_g3a_r_point_after_u, t_cne,
       s_g3a_r_point_cne, 0, 0,
       0, &DataParser::add_text, &DataParser::g3a_r_point_cne);

  init(s_g3a_r_point_after_u, t_cnu,
       s_g3a_r_point_cnu, 0, 0,
       0, &DataParser::add_text, &DataParser::g3a_r_point_cnu);

  init(s_g3a_r_point_after_u, t_cee,
       s_g3a_r_point_cee, 0, 0,
       0, &DataParser::add_text, &DataParser::g3a_r_point_cee);

  init(s_g3a_r_point_after_u, t_ceu,
       s_g3a_r_point_ceu, 0, 0,
       0, &DataParser::add_text, &DataParser::g3a_r_point_ceu);

  init(s_g3a_r_point_after_u, t_cuu,
       s_g3a_r_point_cuu, 0, 0,
       0, &DataParser::add_text, &DataParser::g3a_r_point_cuu);

 // ..... <point> <x/y/z-given/correction/adjusted> ...............

  init(s_g3a_r_point_after_u, t_x_given,
       s_g3a_r_point_x_given, 0, 0,
       0, &DataParser::add_text, &DataParser::g3a_r_point_x_given);

  init(s_g3a_r_point_after_u, t_x_correction,
       s_g3a_r_point_x_correction, 0, 0,
       0, &DataParser::add_text, &DataParser::g3a_r_point_x_correction);

  init(s_g3a_r_point_after_u, t_x_adjusted,
       s_g3a_r_point_x_adjusted, 0, 0,
       0, &DataParser::add_text, &DataParser::g3a_r_point_x_adjusted);

  init(s_g3a_r_point_after_u, t_y_given,
       s_g3a_r_point_y_given, 0, 0,
       0, &DataParser::add_text, &DataParser::g3a_r_point_y_given);

  init(s_g3a_r_point_after_u, t_y_correction,
       s_g3a_r_point_y_correction, 0, 0,
       0, &DataParser::add_text, &DataParser::g3a_r_point_y_correction);

  init(s_g3a_r_point_after_u, t_y_adjusted,
       s_g3a_r_point_y_adjusted, 0, 0,
       0, &DataParser::add_text, &DataParser::g3a_r_point_y_adjusted);

  init(s_g3a_r_point_after_u, t_z_given,
       s_g3a_r_point_z_given, 0, 0,
       0, &DataParser::add_text, &DataParser::g3a_r_point_z_given);

  init(s_g3a_r_point_after_u, t_z_correction,
       s_g3a_r_point_z_correction, 0, 0,
       0, &DataParser::add_text, &DataParser::g3a_r_point_z_correction);

  init(s_g3a_r_point_after_u, t_z_adjusted,
       s_g3a_r_point_z_adjusted, 0, 0,
       0, &DataParser::add_text, &DataParser::g3a_r_point_z_adjusted);

  // ..... <point> <cxx>, <cxy>, <cxz>, <cyy>, <cyz>, <czz> ........

  init(s_g3a_r_point_after_u, t_cxx,
       s_g3a_r_point_cxx, 0, 0,
       0, &DataParser::add_text, &DataParser::g3a_r_point_cxx);

  init(s_g3a_r_point_after_u, t_cxy,
       s_g3a_r_point_cxy, 0, 0,
       0, &DataParser::add_text, &DataParser::g3a_r_point_cxy);

  init(s_g3a_r_point_after_u, t_cxz,
       s_g3a_r_point_cxz, 0, 0,
       0, &DataParser::add_text, &DataParser::g3a_r_point_cxz);

  init(s_g3a_r_point_after_u, t_cyy,
       s_g3a_r_point_cyy, 0, 0,
       0, &DataParser::add_text, &DataParser::g3a_r_point_cyy);

  init(s_g3a_r_point_after_u, t_cyz,
       s_g3a_r_point_cyz, 0, 0,
       0, &DataParser::add_text, &DataParser::g3a_r_point_cyz);

  init(s_g3a_r_point_after_u, t_czz,
       s_g3a_r_point_czz, 0, 0,
       0, &DataParser::add_text, &DataParser::g3a_r_point_czz);

 // ..... <point> <b/l/h-given/correction/adjusted> ...............

  init(s_g3a_r_point_after_u, t_b_given,
       s_g3a_r_point_b_given, 0, 0,
       0, &DataParser::add_text, &DataParser::g3a_r_point_b_given);

  init(s_g3a_r_point_after_u, t_b_correction,
       s_g3a_r_point_b_correction, 0, 0,
       0, &DataParser::add_text, &DataParser::g3a_r_point_b_correction);

  init(s_g3a_r_point_after_u, t_b_adjusted,
       s_g3a_r_point_b_adjusted, 0, 0,
       0, &DataParser::add_text, &DataParser::g3a_r_point_b_adjusted);

  init(s_g3a_r_point_after_u, t_l_given,
       s_g3a_r_point_l_given, 0, 0,
       0, &DataParser::add_text, &DataParser::g3a_r_point_l_given);

  init(s_g3a_r_point_after_u, t_l_correction,
       s_g3a_r_point_l_correction, 0, 0,
       0, &DataParser::add_text, &DataParser::g3a_r_point_l_correction);

  init(s_g3a_r_point_after_u, t_l_adjusted,
       s_g3a_r_point_l_adjusted, 0, 0,
       0, &DataParser::add_text, &DataParser::g3a_r_point_l_adjusted);

  init(s_g3a_r_point_after_u, t_h_given,
       s_g3a_r_point_h_given, 0, 0,
       0, &DataParser::add_text, &DataParser::g3a_r_point_h_given);

  init(s_g3a_r_point_after_u, t_h_correction,
       s_g3a_r_point_h_correction, 0, 0,
       0, &DataParser::add_text, &DataParser::g3a_r_point_h_correction);

  init(s_g3a_r_point_after_u, t_h_adjusted,
       s_g3a_r_point_h_adjusted, 0, 0,
       0, &DataParser::add_text, &DataParser::g3a_r_point_h_adjusted); 

 // ..... <adjusted-observations> .................................

  init(s_g3a_r_adj_results_end, t_adj_observations,
       s_g3a_o_observations, 0, s_g3a_o_observations_end,
       0, 0, 0); 

}


// callback functions for <g3-adjustment-results>

int DataParser::g3a_s_adj_results(const char *name, const char **atts)
{
  no_attributes( name, atts );
  state = next[state][tag(name)];

  delete g3adj->adj;  
  g3adj->adj = new g3::AdjustmentResults;
  
  return 0;
}


int DataParser::g3a_s_adj_results(const char *name)
{
  objects.push_back( new DataObject::g3_adj_results(g3adj->adj) );
  g3adj->adj = 0;

  return  end_tag(name);
}


void DataParser::g3a_text_string(std::string& str)
{
  std::string::const_iterator b=text_buffer.begin();
  std::string::const_iterator e=text_buffer.end();
  TrimWhiteSpaces(b, e);
  str = std::string(b, e);
  text_buffer.clear();
}


void DataParser::g3a_text_float(std::string& str)
{
  std::string::const_iterator b=text_buffer.begin();
  std::string::const_iterator e=text_buffer.end();
  TrimWhiteSpaces(b, e);
  str = std::string(b, e);
  if (!IsFloat(b, e)) error("### bad format of float");
  text_buffer.clear();
}


void DataParser::g3a_text_integer(std::string& str)
{
  std::string::const_iterator b=text_buffer.begin();
  std::string::const_iterator e=text_buffer.end();
  TrimWhiteSpaces(b, e);
  str = std::string(b, e);
  if (!IsInteger(b, e)) error("### bad format of integer");
  text_buffer.clear();
}


int DataParser::g3a_s_algorithm(const char *name)
{
  g3a_text_string(g3adj->adj->algorithm);

  return  end_tag(name);
}


int DataParser::g3a_s_ell_caption(const char *name)
{
  g3a_text_string(g3adj->adj->ell_cap);

  return  end_tag(name);
}


int DataParser::g3a_s_ell_id(const char *name)
{
  g3a_text_string(g3adj->adj->ell_id);
  return  end_tag(name);
}


int DataParser::g3a_s_ell_a(const char *name)
{
  g3a_text_float(g3adj->adj->ell_a);
  return  end_tag(name);
}


int DataParser::g3a_s_ell_b(const char *name)
{
  g3a_text_float(g3adj->adj->ell_b);
  return  end_tag(name);
}


int DataParser::g3a_s_parameters(const char *name)
{
  g3a_text_integer(g3adj->adj->parameters);
  return  end_tag(name);
}


int DataParser::g3a_s_equations(const char *name)
{
  g3a_text_integer(g3adj->adj->equations);
  return  end_tag(name);
}


int DataParser::g3a_s_defect(const char *name)
{
  g3a_text_integer(g3adj->adj->defect);
  return  end_tag(name);
}


int DataParser::g3a_s_redundancy(const char *name)
{
  g3a_text_integer(g3adj->adj->redundancy);
  return  end_tag(name);
}


int DataParser::g3a_s_sum_of_squares(const char *name)
{
  g3a_text_float(g3adj->adj->sum_of_squares);
  return  end_tag(name);
}


int DataParser::g3a_s_apriori_var(const char *name)
{
  g3a_text_float(g3adj->adj->apriori_var);
  return  end_tag(name);
}


int DataParser::g3a_s_aposteriori_var(const char *name)
{
  g3a_text_float(g3adj->adj->aposteriori_var);
  return  end_tag(name);
}


int DataParser::g3a_s_variance_factor(const char *name)
{
  g3a_text_string(g3adj->adj->variance_factor);
  return  end_tag(name);
}


int DataParser::g3a_s_design_m_graph(const char *name)
{
  g3a_text_string(g3adj->adj->design_m_graph);
  return  end_tag(name);
}


int DataParser::g3a_r_point(const char *name, const char **atts)
{
  no_attributes( name, atts );
  state = next[state][tag(name)];

  g3adj->adj->point.clear();
  g3adj->adj->point.n = "unused";
  g3adj->adj->point.e = "unused";
  g3adj->adj->point.u = "unused";
  
  return 0;
}


int DataParser::g3a_r_point(const char *name)
{
  g3adj->adj->points.push_back(g3adj->adj->point);

  return  end_tag(name);
}


int DataParser::g3a_r_point_id(const char *name)
{
  g3a_text_string(g3adj->adj->point.id);
  return  end_tag(name);
}


int DataParser::g3a_r_point_n_type(const char *name)
{
  g3adj->adj->point.n = name;
  return  end_tag(name);
}


int DataParser::g3a_r_point_e_type(const char *name)
{
  g3adj->adj->point.e = name;
  return  end_tag(name);
}


int DataParser::g3a_r_point_u_type(const char *name)
{
  g3adj->adj->point.u = name;
  return  end_tag(name);
}


int DataParser::g3a_r_point_n_dn(const char *name)
{
  g3a_text_float(g3adj->adj->point.n_dn);
  return  end_tag(name);
}


int DataParser::g3a_r_point_e_de(const char *name)
{
  g3a_text_float(g3adj->adj->point.e_de);
  return  end_tag(name);
}


int DataParser::g3a_r_point_u_du(const char *name)
{
  g3a_text_float(g3adj->adj->point.u_du);
  return  end_tag(name);
}


int DataParser::g3a_r_point_n_ind(const char *name)
{
  g3a_text_integer(g3adj->adj->point.n_ind);
  return  end_tag(name);
}


int DataParser::g3a_r_point_e_ind(const char *name)
{
  g3a_text_integer(g3adj->adj->point.e_ind);
  return  end_tag(name);
}


int DataParser::g3a_r_point_u_ind(const char *name)
{
  g3a_text_integer(g3adj->adj->point.u_ind);
  return  end_tag(name);
}


int DataParser::g3a_r_point_cnn(const char *name)
{
  g3a_text_float(g3adj->adj->point.cnn);
  return  end_tag(name);
}


int DataParser::g3a_r_point_cne(const char *name)
{
  g3a_text_float(g3adj->adj->point.cne);
  return  end_tag(name);
}


int DataParser::g3a_r_point_cnu(const char *name)
{
  g3a_text_float(g3adj->adj->point.cnu);
  return  end_tag(name);
}


int DataParser::g3a_r_point_cee(const char *name)
{
  g3a_text_float(g3adj->adj->point.cee);
  return  end_tag(name);
}


int DataParser::g3a_r_point_ceu(const char *name)
{
  g3a_text_float(g3adj->adj->point.ceu);
  return  end_tag(name);
}


int DataParser::g3a_r_point_cuu(const char *name)
{
  g3a_text_float(g3adj->adj->point.cuu);
  return  end_tag(name);
}


int DataParser::g3a_r_point_x_given(const char *name)
{
  g3a_text_float(g3adj->adj->point.x_given);
  return  end_tag(name);
}


int DataParser::g3a_r_point_x_correction(const char *name)
{
  g3a_text_float(g3adj->adj->point.x_correction);
  return  end_tag(name);
}


int DataParser::g3a_r_point_x_adjusted(const char *name)
{
  g3a_text_float(g3adj->adj->point.x_adjusted);
  return  end_tag(name);
}


int DataParser::g3a_r_point_y_given(const char *name)
{
  g3a_text_float(g3adj->adj->point.y_given);
  return  end_tag(name);
}


int DataParser::g3a_r_point_y_correction(const char *name)
{
  g3a_text_float(g3adj->adj->point.y_correction);
  return  end_tag(name);
}


int DataParser::g3a_r_point_y_adjusted(const char *name)
{
  g3a_text_float(g3adj->adj->point.y_adjusted);
  return  end_tag(name);
}


int DataParser::g3a_r_point_z_given(const char *name)
{
  g3a_text_float(g3adj->adj->point.z_given);
  return  end_tag(name);
}


int DataParser::g3a_r_point_z_correction(const char *name)
{
  g3a_text_float(g3adj->adj->point.z_correction);
  return  end_tag(name);
}


int DataParser::g3a_r_point_z_adjusted(const char *name)
{
  g3a_text_float(g3adj->adj->point.x_adjusted);
  return  end_tag(name);
}


int DataParser::g3a_r_point_cxx(const char *name)
{
  g3a_text_float(g3adj->adj->point.cxx);
  return  end_tag(name);
}


int DataParser::g3a_r_point_cxy(const char *name)
{
  g3a_text_float(g3adj->adj->point.cxy);
  return  end_tag(name);
}


int DataParser::g3a_r_point_cxz(const char *name)
{
  g3a_text_float(g3adj->adj->point.cxz);
  return  end_tag(name);
}


int DataParser::g3a_r_point_cyy(const char *name)
{
  g3a_text_float(g3adj->adj->point.cyy);
  return  end_tag(name);
}


int DataParser::g3a_r_point_cyz(const char *name)
{
  g3a_text_float(g3adj->adj->point.cyz);
  return  end_tag(name);
}


int DataParser::g3a_r_point_czz(const char *name)
{
  g3a_text_float(g3adj->adj->point.czz);
  return  end_tag(name);
}


int DataParser::g3a_r_point_b_given(const char *name)
{
  g3a_text_string(g3adj->adj->point.b_given);
  return  end_tag(name);
}


int DataParser::g3a_r_point_b_correction(const char *name)
{
  g3a_text_string(g3adj->adj->point.b_correction);
  return  end_tag(name);
}


int DataParser::g3a_r_point_b_adjusted(const char *name)
{
  g3a_text_string(g3adj->adj->point.b_adjusted);
  return  end_tag(name);
}


int DataParser::g3a_r_point_l_given(const char *name)
{
  g3a_text_string(g3adj->adj->point.l_given);
  return  end_tag(name);
}


int DataParser::g3a_r_point_l_correction(const char *name)
{
  g3a_text_string(g3adj->adj->point.l_correction);
  return  end_tag(name);
}


int DataParser::g3a_r_point_l_adjusted(const char *name)
{
  g3a_text_string(g3adj->adj->point.l_adjusted);
  return  end_tag(name);
}


int DataParser::g3a_r_point_h_given(const char *name)
{
  g3a_text_string(g3adj->adj->point.h_given);
  return  end_tag(name);
}


int DataParser::g3a_r_point_h_correction(const char *name)
{
  g3a_text_string(g3adj->adj->point.h_correction);
  return  end_tag(name);
}


int DataParser::g3a_r_point_h_adjusted(const char *name)
{
  g3a_text_string(g3adj->adj->point.h_adjusted);
  return  end_tag(name);
}