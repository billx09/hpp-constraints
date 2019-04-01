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

#ifndef HPP_CONSTRAINTS_STATIC_TORQUE_HH
# define HPP_CONSTRAINTS_STATIC_TORQUE_HH

# include <hpp/constraints/differentiable-function.hh>
# include <hpp/constraints/matrix-view.hh>

namespace hpp {
  namespace constraints {

    /// \addtogroup constraints
    /// \{

    HPP_PREDEF_CLASS (StaticTorque);
    typedef boost::shared_ptr<StaticTorque> StaticTorquePtr_t;

    /// Compute the static torques of the robot.
    /// \todo add an API to provide external forces.
    class HPP_CONSTRAINTS_DLLAPI StaticTorque : public DifferentiableFunction
    {
    public:
      virtual ~StaticTorque () {}

      static StaticTorquePtr_t create (const DevicePtr_t& device, const std::string& name);

      void selection (const segments_t& idxs);

      /// Whether to return the absolute value minus the torque limit.
      /// If set to true, then checking the bound limits breaks down to
      /// \f$ StaticTorque(q) \le 0 \f$
      void absoluteValue (bool set);

      /// Display object in a stream
      virtual std::ostream& print (std::ostream& o) const;

    protected:
      /// \param name function's name
      StaticTorque (const DevicePtr_t& device, const std::string& name);

      /// User implementation of function evaluation
      virtual void impl_compute (LiegroupElementRef result, vectorIn_t argument) const;

      virtual void impl_jacobian (matrixOut_t jacobian, vectorIn_t arg) const;

    private:
      DevicePtr_t device_;
      /// The selection of DoFs
      Eigen::RowBlockIndices jView_;
      /// Whether to return the absolute value minus the torque limit
      bool absVal_;
    }; // class StaticTorque

    ImplicitPtr_t torqueLimits (const DevicePtr_t& device, const std::string& name);

    /// \}
  } // namespace constraints
} // namespace hpp


#endif // HPP_CONSTRAINTS_STATIC_TORQUE_HH
