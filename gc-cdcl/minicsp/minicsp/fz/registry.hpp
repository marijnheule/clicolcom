/* -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */
/*****************************************************************************
  minicsp

Copyright 2010 George Katsirelos

Minicsp is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

Minicsp is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with minicsp.  If not, see <http://www.gnu.org/licenses/>.

Most of this file was distributed under the following license as part
of Gecode.

  Main authors:
     Guido Tack <tack@gecode.org>

  Copyright:
     Guido Tack, 2007

  Last modified:
     $Date: 2006-12-11 03:27:31 +1100 (Mon, 11 Dec 2006) $ by $Author: schulte $
     $Revision: 4024 $

  This file is part of Gecode, the generic constraint
  development environment:
     http://www.gecode.org

  Permission is hereby granted, free of charge, to any person obtaining
  a copy of this software and associated documentation files (the
  "Software"), to deal in the Software without restriction, including
  without limitation the rights to use, copy, modify, merge, publish,
  distribute, sublicense, and/or sell copies of the Software, and to
  permit persons to whom the Software is furnished to do so, subject to
  the following conditions:

  The above copyright notice and this permission notice shall be
  included in all copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
  OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*****************************************************************************/

#ifndef __MINICSP_FLATZINC_REGISTRY_HH__
#define __MINICSP_FLATZINC_REGISTRY_HH__

#include "flatzinc.hpp"
#include <string>
#include <map>

namespace minicsp {

namespace FlatZinc {

  /// Map from constraint identifier to constraint posting functions
  class Registry {
  public:
    /// Type of constraint posting function
    typedef void (*poster) (Solver&,
                            FlatZincModel&,
                            const ConExpr&,
                            AST::Node*);
    /// Add posting function \a p with identifier \a id
    void add(const std::string& id, poster p);
    /// Post constraint specified by \a ce
    void post(Solver& s, FlatZincModel &m,
              const ConExpr& ce, AST::Node* ann);

  private:
    /// The actual registry
    std::map<std::string,poster> r;
  };

  /// Return global registry object
  Registry& registry(void);

}

}

#endif

// STATISTICS: flatzinc-any
