/*  
    GNU Gama -- adjustment of geodetic networks
    Copyright (C) 2002  Ales Cepek <cepek@fsv.cvut.cz>

    This file is part of the GNU Gama C++ library.
    
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
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/*
 *  $Id: intlist.h,v 1.3 2003/03/25 12:38:33 cepek Exp $
 */

#ifndef GNU_gama___GaMaLib_Integer_list____GaMaLib_Integer_list
#define GNU_gama___GaMaLib_Integer_list____GaMaLib_Integer_list

#include <cstddef>

namespace GNU_gama {

template <class Index=std::size_t> 

  class IntegerList {     // Integer List class
    
    IntegerList (const IntegerList&); 
    void operator=(const IntegerList&);

    Index*  m;
    Index*  e;

    public:
    
    IntegerList() : m(0), e(0)
    {
    }
    IntegerList(Index n)
    {
      m = new Index[n];
      e = m + n;
    }
    
    ~IntegerList() 
    {
      delete[] m;
    }

    Index dim() const { return Index(e-m); }

    void reset()
    {
      delete[] m;
      m = e = 0;
    }
    void reset(Index n)
    {
      delete[] m;
      m = new Index[n];
      e = m + n;
    }

    typedef Index* iterator;
    typedef const Index* const_iterator;
    iterator begin() { return m; }
    iterator end()   { return e; } 
    const_iterator begin() const { return m; }
    const_iterator end()   const { return e; } 
    
  };
 
}   // namespace GNU_gama

#endif



#ifdef GNU_gama_sparse_demo

#include <iostream>

using namespace std;
using namespace GNU_gama;


int main()
{
  cout << "\n---  Integer List demo  ------------------------------------\n\n";

  IntegerList<> ilist;
  ilist.reset(20);
  IntegerList<>::iterator c = ilist.begin();
  IntegerList<>::const_iterator e = ilist.end();
  for (int i=1, j=1, k=1; c!=e; k = i+j, i=j, j=k)
    {
      *c++ = i;
    }

  for (IntegerList<>::const_iterator a=ilist.begin(), b=ilist.end(); a!=b; a++)
    {
      cout << *a << " ";
    }
  cout << "\n\n";
}


#endif