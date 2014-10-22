// Copyright (c) 2014, LAAS-CNRS
// Authors: Joseph Mirabel (joseph.mirabel@laas.fr)
//
// This file is part of hpp-constraints.
// hpp-constraints is free software: you can redistribute it
// and/or modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation, either version
// 3 of the License, or (at your option) any later version.
//
// hpp-constraints is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Lesser Public License for more details.  You should have
// received a copy of the GNU Lesser General Public License along with
// hpp-constraints. If not, see <http://www.gnu.org/licenses/>.

#ifndef HPP_CONSTRAINTS_STATIC_STABILITY_HH
# define HPP_CONSTRAINTS_STATIC_STABILITY_HH

# include <vector>
# include <hpp/core/differentiable-function.hh>

# include "hpp/constraints/fwd.hh"

namespace hpp {
  namespace constraints {
    class Triangle {
      public:
        /// Represent a triangle.
        Triangle (const fcl::Vec3f& p0, const fcl::Vec3f& p1, const fcl::Vec3f& p2):
          p0_ (p0), p1_ (p1), p2_ (p2)
        {
          n_ = (p1_ - p0_).cross (p2_ - p0_);  
          assert (!n_.isZero ());
          n_.normalize ();
          c_ = ( p0_ + p1_ + p2_ ) / 3;
          n0_ = (p2_ - p1_).cross (n_); n0_.normalize ();
          n1_ = (p0_ - p2_).cross (n_); n0_.normalize ();
          n2_ = (p1_ - p0_).cross (n_); n0_.normalize ();
          nxn0_ = n_.cross (n0_);
        }

        /// Intersection with a line defined by a point and a vector.
        inline fcl::Vec3f intersection (const fcl::Vec3f& A, const fcl::Vec3f& u) const {
          assert (std::abs (n_.dot (u)) > 1e-8);
          return A + u * (n_.dot(c_ - A)) / n_.dot (u);
        }

        /// Check whether the intersection of the line defined by A and u
        /// onto the plane containing the triangle is inside the triangle.
        inline bool isInside (const fcl::Vec3f& A, const fcl::Vec3f& u) const {
          fcl::Vec3f Ap = intersection (A, u);
          return isInside (Ap);
        }
        inline bool isInside (const fcl::Vec3f& Ap) const {
          if (n2_.dot (Ap-p0_) >= 0) return false;
          if (n0_.dot (Ap-p1_) >= 0) return false;
          if (n1_.dot (Ap-p2_) >= 0) return false;
          return true;
        }

        /// Return the shortest distance from a point to the triangle.
        /// A negative value means the point is inside the triangle.
        /// \param A a point already in the plane containing the triangle.
        inline value_type distance (const fcl::Vec3f& A) {
          value_type d0 = dist (A-p0_, p1_-p0_, n2_),
                     d1 = dist (A-p1_, p2_-p1_, n0_),
                     d2 = dist (A-p2_, p0_-p2_, n1_);
          /// We should return the smallest in norm.
          if (std::abs (d0) > std::abs (d1)) {
            if (std::abs (d1) > std::abs (d2)) return d2;
            else return d1;
          }
          if (std::abs (d0) > std::abs (d2)) return d2;
          else return d0;
        }

        inline const fcl::Vec3f& planeXaxis () const { return   n0_; }
        inline const fcl::Vec3f& planeYaxis () const { return nxn0_; }
        inline const fcl::Vec3f& normal () const { return n_; }
        inline const fcl::Vec3f& center () const { return c_; }

      private:
        /// Return the distance between the point A and the segment
        /// [P, v] oriented by u.
        /// w = PA.
        inline value_type dist (const fcl::Vec3f& w, const fcl::Vec3f& v, const fcl::Vec3f& u) const {
          value_type c1, c2;
          c1 = v.dot (w);
          if (c1 <= 0)
            return (u.dot (w) > 0)?(w.norm()):(- w.norm());
          c2 = v.norm ();
          if (c2 <= c1)
            return (u.dot (w) > 0)?((w - v).norm()):((v-w).norm());
          return u.dot (w);
        }

        fcl::Vec3f p0_, p1_, p2_, n_, c_;
        /// n_i is the vector of norm 1 perpendicular to
        /// P_{i+1}P_i and n_.
        fcl::Vec3f n0_, n1_, n2_, nxn0_;
    };

    class StaticStabilityGravity : public DifferentiableFunction {
      public:
        static fcl::Vec3f gravity;

        /// Constructor
        /// \param robot the robot the constraints is applied to,
        /// \param joint the joint to which the object is attached,
        /// \param com COM of the object in the joint frame.
        StaticStabilityGravity (const DevicePtr_t& robot, const JointPtr_t& joint,
            const fcl::Vec3f& com);

        static StaticStabilityGravityPtr_t create (const DevicePtr_t& robot,
            const JointPtr_t& joint, const fcl::Vec3f& com);

        void addObjectTriangle (const Triangle& t);

        void addFloorTriangle (const Triangle& t);

      private:
        void impl_compute (vectorOut_t result, ConfigurationIn_t argument) const;

        void impl_jacobian (matrixOut_t jacobian, ConfigurationIn_t argument) const;

        void selectTriangles () const;

        DevicePtr_t robot_;
        JointPtr_t joint_;
        fcl::Vec3f com_;
        mutable fcl::Vec3f comGlobalFrame_;

        typedef std::vector <Triangle> Triangles;
        /// Triangles with coordinates expressed in joint frame.
        Triangles objectTriangles_;
        mutable Triangles::const_iterator object_;
        /// Triangles with coordinates expressed in world frame.
        Triangles floorTriangles_;
        mutable Triangles::const_iterator floor_;

        mutable Eigen::Matrix <double, 3, 3> Rcx_;
        mutable vector_t n_;
    };
  } // namespace constraints
} // namespace hpp

#endif //  HPP_CONSTRAINTS_STATIC_STABILITY_HH