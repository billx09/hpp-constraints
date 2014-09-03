//
// Copyright (c) 2014 CNRS
// Authors: Florent Lamiraux
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

#include <hpp/model/device.hh>
#include <hpp/model/fcl-to-eigen.hh>
#include <hpp/model/joint.hh>
#include <hpp/constraints/relative-position.hh>

namespace hpp {
  namespace constraints {
    static size_type size (std::vector<bool> mask)
    {
      size_type res = 0;
      for (std::vector<bool>::iterator it = mask.begin (); it != mask.end ();
	   ++it) {
	if (*it) ++res;
      }
      return res;
    }


    RelativePositionPtr_t RelativePosition::create
    (const DevicePtr_t& robot, const JointPtr_t&  joint1,
     const JointPtr_t& joint2, const vector3_t& point1,
     const vector3_t& point2, std::vector <bool> mask)
    {
      RelativePosition* ptr = new RelativePosition
	(robot, joint1, joint2, point1, point2, mask);
      RelativePositionPtr_t shPtr (ptr);
      return shPtr;
    }

    static void cross (const vector3_t& v, eigen::matrix3_t& m)
    {
      m (0,1) = -v [2]; m (1,0) = v [2];
      m (0,2) = v [1]; m (2,0) = -v [1];
      m (1,2) = -v [0]; m (2,1) = v [0];
    }

    RelativePosition::RelativePosition (const DevicePtr_t& robot,
					const JointPtr_t& joint1,
					const JointPtr_t& joint2,
					const vector3_t& point1,
					const vector3_t& point2,
					std::vector <bool> mask) :
      DifferentiableFunction (robot->configSize (), robot->numberDof (),
		size (mask), "RelativePosition"),
      robot_ (robot), joint1_ (joint1), joint2_ (joint2),
      point1_ (point1), point2_ (point2), mask_ (mask),
      jacobian_ (3, robot->numberDof ())
    {
      cross1_.setZero ();
      cross2_.setZero ();
      result_.resize (3); result_.setZero ();
      jacobian_.setZero ();
    }

     void RelativePosition::impl_compute
     (vectorOut_t result, ConfigurationIn_t argument) const throw ()
    {
      robot_->currentConfiguration (argument);
      robot_->computeForwardKinematics ();
      const Transform3f& M1 = joint1_->currentTransformation ();
      const Transform3f& M2 = joint2_->currentTransformation ();
      global1_ = M1.transform (point1_);
      global2_ = M2.transform (point2_);
      model::toEigen (global1_ - global2_, result_);
      size_type index = 0;
      if (mask_ [0]) {
	result [index] = result_ [0]; ++index;
      }
      if (mask_ [1]) {
	result [index] = result_ [1]; ++index;
      }
      if (mask_ [2]) {
	result [index] = result_ [2]; ++index;
      }
    }

  void RelativePosition::impl_jacobian (matrixOut_t jacobian,
					ConfigurationIn_t arg) const throw ()
    {
      robot_->currentConfiguration (arg);
      robot_->computeForwardKinematics ();
      const Transform3f& M1 = joint1_->currentTransformation ();
      const Transform3f& M2 = joint2_->currentTransformation ();
      const JointJacobian_t& Jjoint1 (joint1_->jacobian ());
      const JointJacobian_t& Jjoint2 (joint2_->jacobian ());
      R1x1_ = M1.getRotation () * point1_;
      // -[R1 (q) x1]x
      cross (-R1x1_, cross1_);
      R2x2_ = M2.getRotation () * point2_;
      // [R2 (q) x2]x
      cross (R2x2_, cross2_);
      jacobian_.leftCols (Jjoint1.cols ()) = cross1_ * Jjoint1.bottomRows (3)
	+ Jjoint1.topRows (3) + cross2_ * Jjoint2.bottomRows (3)
	- Jjoint2.topRows (3);
      size_type index = 0;
      if (mask_ [0]) {
	jacobian.row (index) = jacobian_.row (0); ++index;
      }
      if (mask_ [1]) {
	jacobian.row (index) = jacobian_.row (1); ++index;
      }
      if (mask_ [2]) {
	jacobian.row (index) = jacobian_.row (2); ++index;
      }
    }

  } // namespace constraints
} // namespace hpp
