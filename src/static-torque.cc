//
// Copyright (c) 2019 CNRS
// Authors: Joseph Mirabel
//
// This file is part of hpp-constraints
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
// hpp-constraints  If not, see
// <http://www.gnu.org/licenses/>.

#include <hpp/constraints/static-torque.hh>

#include <pinocchio/algorithm/rnea.hpp>
#include <pinocchio/algorithm/rnea-derivatives.hpp>

#include <hpp/pinocchio/device.hh>
#include <hpp/pinocchio/joint-collection.hh>

#include <hpp/constraints/implicit.hh>

namespace hpp {
  namespace constraints {
    StaticTorquePtr_t StaticTorque::create (const DevicePtr_t& device, const std::string& name)
    {
      return StaticTorquePtr_t (new StaticTorque (device, name));
    }

    void StaticTorque::selection (const segments_t& idxs)
    {
      jView_ = Eigen::RowBlockIndices(idxs);
      outputSpace_ = LiegroupSpace::Rn(jView_.nbRows());
      //TODO update the active parameters.
    }

    void StaticTorque::absoluteValue (bool set)
    {
      absVal_ = set;
    }

    std::ostream& StaticTorque::print (std::ostream& o) const
    {
      return o << "StaticTorque: " << name () << " "
        << (absVal_ ? "absolute value" : "") << incindent << iendl
        << "selection: " << jView_
        << decindent;
    }

    StaticTorque::StaticTorque (const DevicePtr_t& device, const std::string& name)
      : DifferentiableFunction (device->configSize(), device->numberDof(),
          device->numberDof(), name),
      device_ (device)
    {
    }

    void StaticTorque::impl_compute (LiegroupElementRef result, vectorIn_t q) const
    {
      hpp::pinocchio::DeviceSync device (device_);
      const hpp::pinocchio::Model& model = device.model();
      hpp::pinocchio::Data& data = device.data();
      // If velocity and accelerations are not zero, one should use RNEA instead.
      ::pinocchio::computeGeneralizedGravity(model, data, q);
      if (absVal_)
        result.vector().noalias() = jView_.rview(data.g.cwiseAbs() - model.effortLimit);
      else
        result.vector() = jView_.rview (data.g);
    }

    void StaticTorque::impl_jacobian (matrixOut_t jacobian, vectorIn_t q) const
    {
      hpp::pinocchio::DeviceSync device (device_);
      const hpp::pinocchio::Model& model = device.model();
      hpp::pinocchio::Data& data = device.data();
      // If velocity and accelerations are not zero, one should use RNEA instead.
      data.dtau_dq.setZero();
      ::pinocchio::computeGeneralizedGravityDerivatives(model, data, q, data.dtau_dq);

      if (absVal_) {
        ::pinocchio::computeGeneralizedGravity(device.model(), data, q);
        jacobian = jView_.rview((data.g.array() > 0).select (data.dtau_dq, -data.dtau_dq));
      } else
        jacobian = jView_.rview(data.dtau_dq);
    }

    ImplicitPtr_t torqueLimits (const DevicePtr_t& device, const std::string& name)
    {
      const pinocchio::Model& model = device->model();

      // Compute the DoFs
      ArrayXb torqueBoundIsNotNaN = (model.effortLimit.array() == model.effortLimit.array());
      ArrayXb boundedTorque = torqueBoundIsNotNaN
        && (model.effortLimit.array() < std::numeric_limits<value_type>::max());

      StaticTorquePtr_t function = StaticTorque::create (device, name);
      function->selection (Eigen::BlockIndex::fromLogicalExpression(boundedTorque));
      function->absoluteValue (true);

      ImplicitPtr_t constraint = Implicit::create (function,
          ComparisonTypes_t (function->outputSize(), Inferior));
      return constraint;
    }

  } // namespace constraints
} // namespace hpp
