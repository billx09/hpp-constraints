// Copyright (c) 2019, LAAS-CNRS
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

#include <hpp/constraints/static-torque.hh>

#define BOOST_TEST_MODULE static_torque
#include <boost/test/included/unit_test.hpp>

#include <pinocchio/multibody/model.hpp>
#include <hpp/pinocchio/liegroup-space.hh>
#include <hpp/pinocchio/simple-device.hh>
#include <hpp/constraints/implicit.hh>

BOOST_AUTO_TEST_CASE(static_torque)
{
  using namespace hpp::pinocchio;
  using namespace hpp::constraints;
  DevicePtr_t device = unittest::makeDevice (
        //unittest::ManipulatorArm2);
        unittest::HumanoidRomeo);
  ImplicitPtr_t tauLimits = torqueLimits (device, "torque_limit");

  LiegroupSpacePtr_t expectedSpace = LiegroupSpace::Rn(device->numberDof() - 6);
  BOOST_CHECK_MESSAGE (*tauLimits->function().outputSpace() == *expectedSpace,
      "Function\n" << tauLimits->function() << "\nThe output space should be "
      << *expectedSpace << " instead of " << *tauLimits->function().outputSpace()
      << ".\nThe effort limits " << device->model().effortLimit.transpose()
      << " should be invalid (NaN, Inf or " << std::numeric_limits<value_type>::max()
      << ") on the first 6 dofs and finite values for the rest.");
}
