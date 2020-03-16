//
// Copyright (c) 2020 CNRS
// Authors: Joseph Mirabel
//
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

#ifndef HPP_CONSTRAINTS_CONVEX_COLLISION_AVOIDANCE_HH
# define HPP_CONSTRAINTS_CONVEX_COLLISION_AVOIDANCE_HH

# include <hpp/pinocchio/liegroup-element.hh>

# include <hpp/fcl/narrowphase/gjk.h>
# include <hpp/constraints/fwd.hh>
# include <hpp/constraints/differentiable-function.hh>

namespace hpp {
namespace constraints {

class ConvexCollisionAvoidance;
typedef boost::shared_ptr <ConvexCollisionAvoidance> ConvexCollisionAvoidancePtr_t;

/// ConvexCollisionAvoidance sets of objects
///
/// This function maps to a configuration of a robot, the distance
///   \li either between objects of a joints and objects of another joint,
///   \li or objects of a joint with a list of fixed objects.
///
/// The above type of distance is determined by the method "create" called.
class HPP_CONSTRAINTS_DLLAPI ConvexCollisionAvoidance :
  public DifferentiableFunction
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  /// Create instance and return shared pointer
  ///
  /// \param name name of the constraint,
  /// \param robot robot that own the bodies,
  static ConvexCollisionAvoidancePtr_t create (const std::string& name,
                                            const DevicePtr_t& robot);

  virtual ~ConvexCollisionAvoidance () {}

  void addCollisionPair (::pinocchio::GeomIndex geom1,
      ::pinocchio::GeomIndex geom2,
      value_type distanceUpperBound, value_type precision = 1e-3);

  void addCollisionPair (::pinocchio::PairIndex index,
      value_type distanceUpperBound, value_type precision = 1e-3);

  void addCollisionPair (const std::string& geom1, const std::string& geom2,
      value_type distanceUpperBound, value_type precision = 1e-3);

  void addObstacleModel (const pinocchio::GeomModelPtr_t& model,
      const pinocchio::GeomDataPtr_t& data)
  {
    obsModel_ = model;
    obsData_ = data;
  }

protected:
  /// Protected constructor
  ///
  /// \param name name of the constraint,
  /// \param robot robot that own the bodies,
  ConvexCollisionAvoidance (const std::string& name,
      const DevicePtr_t& robot);

  virtual void impl_compute (LiegroupElementRef result,
                             ConfigurationIn_t argument) const;
  virtual void impl_jacobian (matrixOut_t jacobian,
                              ConfigurationIn_t arg) const;
private:
  void addCollisionPair (
      bool robot1, ::pinocchio::GeomIndex geom1,
      bool robot2, ::pinocchio::GeomIndex geom2,
      value_type distanceUpperBound, value_type precision = 1e-3);

  void computeCollision(ConfigurationIn_t arg) const;

  DevicePtr_t robot_;
  pinocchio::GeomModelPtr_t obsModel_;
  pinocchio::GeomDataPtr_t obsData_;

  value_type distanceUpperBound_;

  typedef hpp::fcl::details::MinkowskiDiff MinkowskiDiff;
  typedef hpp::fcl::details::GJK           GJK;
  typedef hpp::fcl::details::EPA           EPA;

  struct Data {
    ::pinocchio::GeomIndex ida, idb;
    bool aRob, bRob;

    MinkowskiDiff shape;
    GJK gjk;
    EPA* epa;

    GJK::Status gjk_status;
    EPA::Status epa_status;
    value_type distanceUpperBound;
    value_type distance;
    /// Closest points and normal are expressed in the frame of the first model.
    vector3_t w0, w1, normal;

    inline Data (value_type precision = 1e-3);
  };
  mutable std::vector<Data> datas_;

  mutable Configuration_t latestArgument_;
  mutable LiegroupElement latestResult_;
}; // class ConvexCollisionAvoidance

} // namespace constraints
} // namespace hpp

#endif // HPP_CONSTRAINTS_CONVEX_COLLISION_AVOIDANCE_HH
