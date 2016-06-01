//
// Copyright (c) 2016 CNRS
// Authors: Joseph Mirabel
//
// This file is part of hpp-constraints.
// hpp-constraints is free software: you can redistribute it
// and/or modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation, either version
// 3 of the License, or (at your option) any later version.
//
// hpp-constraints is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Lesser Public License for more details. You should have
// received a copy of the GNU Lesser General Public License along with
// hpp-constraints. If not, see
// <http://www.gnu.org/licenses/>.

#ifndef HPP_CONSTRAINTS_INTERVALS_HH
# define HPP_CONSTRAINTS_INTERVALS_HH

# include <utility>
# include <hpp/constraints/fwd.hh>
# include <boost/icl/interval_set.hpp>

namespace hpp {
  namespace constraints {
    typedef std::pair<size_type, size_type> SizeInterval_t;
  } // namespace constraints
} // namespace hpp

namespace boost {
  namespace icl {
    template<>
      struct interval_traits< hpp::constraints::SizeInterval_t >
      {
        typedef hpp::constraints::SizeInterval_t interval_type;
        typedef int            domain_type;
        typedef std::less<int> domain_compare;

        static interval_type construct(const domain_type& lo, const domain_type& up) 
        { return interval_type(lo, up - lo); }

        static domain_type lower(const interval_type& inter){ return inter.first; };
        static domain_type upper(const interval_type& inter){ return inter.first + inter.second; };
      };

    template<>
      struct interval_bound_type<hpp::constraints::SizeInterval_t>
      {
        typedef interval_bound_type type;
        BOOST_STATIC_CONSTANT(bound_type, value = interval_bounds::static_right_open);//[lo..up)
      };

  } // namespace icl
} // namespace boost

namespace hpp {
  namespace constraints {
    typedef boost::icl::interval_set<size_type, std::less, SizeInterval_t> SizeIntervals_t;
  } // namespace constraints
} // namespace hpp


#endif // HPP_CONSTRAINTS_INTERVALS_HH
