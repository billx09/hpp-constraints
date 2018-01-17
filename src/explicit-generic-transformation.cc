// Copyright (c) 2016, Joseph Mirabel
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

#include <hpp/constraints/explicit-generic-transformation.hh>

#include <hpp/fcl/math/transform.h>

#include <pinocchio/multibody/model.hpp>

#include <hpp/util/indent.hh>

#include <hpp/pinocchio/configuration.hh>
#include <hpp/pinocchio/device.hh>
#include <hpp/pinocchio/joint.hh>

#include <hpp/constraints/tools.hh>
#include <hpp/constraints/macros.hh>
#include <hpp/constraints/matrix-view.hh>

namespace hpp {
  namespace constraints {
    namespace {
      typedef JointJacobian_t::ConstNRowsBlockXpr<3>::Type HalfJacobian_t;
      inline HalfJacobian_t omega(const JointJacobian_t& j) { return j.bottomRows<3>(); }
      inline HalfJacobian_t trans(const JointJacobian_t& j) { return j.topRows<3>(); }

      template <bool ComputePosition, bool ComputeOrientation> struct output_space;
      template <> struct output_space <true , true > { static inline LiegroupSpacePtr_t get () { return pinocchio::LiegroupSpace::SE3(); } };
      template <> struct output_space <true , false> { static inline LiegroupSpacePtr_t get () { return pinocchio::LiegroupSpace::R3 (); } };
      // template <> struct output_space <false, true > { static inline LiegroupSpacePtr_t get () { return pinocchio::LiegroupSpace::SO3(); } };

      void inputVariable (JointConstPtr_t joint, ArrayXb& conf, ArrayXb& vel)
      {
        while (joint && joint->index() != 0) {
          size_type s = joint->rankInConfiguration(),
                    l = joint->configSize();
          conf.segment(s,l) = conf.segment(s,l).cwiseEqual(false); // Negation

          s = joint->rankInVelocity();
          l = joint->numberDof();
          vel.segment(s,l) = vel.segment(s,l).cwiseEqual(false); // Negation

          hppDout (info, "Adding joint " << joint->name ()
                   << " as input variable.");
          joint = joint->parentJoint();
        }
      }
      template <bool pos, bool ori>
      BlockIndex::segments_t jointConfInterval (JointConstPtr_t j) {
        return BlockIndex::segments_t(1, BlockIndex::segment_t (
              j->rankInConfiguration() + (pos ? 0 : 3),
              j->configSize() - (pos ? 0 : 3) - (ori ? 0 : 4)));
      }
      template <bool pos, bool ori>
      BlockIndex::segments_t jointVelInterval (JointConstPtr_t j) {
        return BlockIndex::segments_t(1, BlockIndex::segment_t (
              j->rankInVelocity() + (pos ? 0 : 3),
              j->numberDof() - (pos ? 0 : 3) - (ori ? 0 : 3)));
      }
    }

    template <int _Options> std::ostream&
      ExplicitGenericTransformation<_Options>::print (std::ostream& os) const
    {
      os << "Explicit" << (IsRelative ? "Relative" : "") <<
            (ComputePosition ? (ComputeOrientation ? "Transformation" : "Position")
                   : "Orientation") << ": " << name()
        << incindent
        << iendl << "Joint1: "          << ((IsRelative && joint1_) ? joint1_->name() : "World")
        << iendl << "Joint2: "          << joint2_->name()
        << iendl << "Reference: "       << incindent
        << iendl << pretty_print (ref_) << decindent
        << iendl << "mask: ";
      for (size_type i=0; i<MaskSize; ++i) os << mask_ [i] << ", ";
      return os << decindent;
    }

    template <int _Options> typename ExplicitGenericTransformation<_Options>::Ptr_t
      ExplicitGenericTransformation<_Options>::create
    (const std::string& name, const DevicePtr_t& robot,
     const JointConstPtr_t& joint1, const JointConstPtr_t& joint2,
     const Transform3f& reference, std::vector <bool> mask)
    {
      ArrayXb conf (ArrayXb::Constant(robot->configSize(), false));
      ArrayXb vel  (ArrayXb::Constant(robot->numberDof(),  false));
      inputVariable (joint1, conf, vel);
      inputVariable (joint2->parentJoint(), conf, vel);

      ExplicitGenericTransformation<_Options>* ptr =
        new ExplicitGenericTransformation<_Options> (name, robot, joint1, joint2,
            reference,
            BlockIndex::fromLogicalExpression (conf),
            jointConfInterval<ComputePosition,ComputeOrientation>(joint2),
            BlockIndex::fromLogicalExpression (vel),
            jointVelInterval<ComputePosition,ComputeOrientation>(joint2),
            mask);
      Ptr_t shPtr (ptr);
      ptr->init (shPtr);
      return shPtr;
    }

    template <int _Options> typename ExplicitGenericTransformation<_Options>::Ptr_t
      ExplicitGenericTransformation<_Options>::create
    (const std::string& name, const DevicePtr_t& robot,
     const JointConstPtr_t& joint1, const JointConstPtr_t& joint2,
     const Transform3f& frame1, const Transform3f& frame2,
     std::vector <bool> mask)
    {
      return create (name, robot, joint1, joint2, frame1 * frame2.inverse(), mask);
    }

    template <int _Options>
      ExplicitGenericTransformation<_Options>::ExplicitGenericTransformation
    (const std::string& name      , const DevicePtr_t& robot,
     const JointConstPtr_t& joint1, const JointConstPtr_t& joint2,
     const Transform3f& ref,
     const segments_t inConf , const segments_t outConf,
     const segments_t inVel  , const segments_t outVel ,
     std::vector <bool> mask):
      DifferentiableFunction (
              BlockIndex::cardinal(inConf),  BlockIndex::cardinal(inVel),
              output_space<ComputePosition, ComputeOrientation>::get (), name),
        robot_ (robot),
        joint1_ (joint1), joint2_ (joint2),
        parentJoint_ (joint2->parentJoint ()),
        inConf_ (inConf),   inVel_  (inVel),
        outConf_ (outConf), outVel_ (outVel),
        mask_ (mask)
    {
      reference (ref);
      assert(IsRelative);
      assert(mask.size()==MaskSize);
      for (int i = 0; i < MaskSize; ++i) { assert (mask[i]); }
    }

    template <int _Options>
    void ExplicitGenericTransformation<_Options>::reference (const Transform3f& reference)
    {
      if (ComputePosition) ref_.translation() = reference.translation();
      else                 ref_.translation().setZero();
      if (ComputeOrientation) ref_.rotation() = reference.rotation();
      else                    ref_.rotation().setIdentity();
    }

    template <int _Options>
    inline void ExplicitGenericTransformation<_Options>::forwardKinematics (const vectorIn_t& arg) const
    {
      qsmall_ = inConf_.rview(robot_->currentConfiguration());
      if (qsmall_ != arg) {
        q_ = robot_->currentConfiguration();
        inConf_.lview(q_) = arg;
        robot_->currentConfiguration(q_);
      }
      robot_->computeForwardKinematics ();
    }

    template <int _Options>
    void ExplicitGenericTransformation<_Options>::impl_compute
    (LiegroupElement& result, ConfigurationIn_t arg) const throw ()
    {
      forwardKinematics (arg);

      // J1 * M1/J1 = J2 * M2/J2
      // J2 = J1 * M1/J1 * M2/J2^{-1}
      // J2 = J2_{parent} * T
      // T = J2_{parent}^{-1} * J2
      // T = J2_{parent}^{-1} * J1 * F1/J1 * F2/J2^{-1}
      freeflyerPose_ = M1() * ref_;

      if (parentJoint_)
        freeflyerPose_ = Mp().actInv(freeflyerPose_);

      freeflyerPose_ = pM2().actInv (freeflyerPose_);

      if (ComputePosition) {
        result.vector ().head<3>() = freeflyerPose_.translation();
        for (size_type i=0; i<3; ++i) { assert(mask_[i]); }
          // if (mask_ [i])
            // result.vector () [index] = freeflyerPose_.translation()[i];
      }
      if (ComputeOrientation) {
        typedef Transform3f::Quaternion_t Q_t;
        result.vector ().tail<4>() = Q_t(freeflyerPose_.rotation()).coeffs();
      }
    }

    template <int _Options>
    void ExplicitGenericTransformation<_Options>::impl_jacobian
    (matrixOut_t jacobian, ConfigurationIn_t arg) const throw ()
    {
      LiegroupElement result (outputSpace ());
      impl_compute (result, arg);

      q_ = robot_->currentConfiguration ();
      outConf_.lview (q_) = result.vector ();
      robot_->currentConfiguration (q_);
      robot_->computeForwardKinematics ();

      cross1_ = se3::skew((R1() * ref_.translation()).eval());
      if (parentJoint_) {
        cross2_.noalias() = se3::skew((pt2() - t1()).eval());
        J2_parent_minus_J1_.noalias() = Jp() - J1();
      } else {
        cross2_.noalias() = - se3::skew(t1());
      }

      // Express velocity of J1 * M1/J1 * M2/J2^{-1} in J2_{parent}.
      if (ComputePosition) {
        if (parentJoint_) {
          tmpJac_.noalias() = (pR2().transpose() * Rp().transpose()) *
            (   cross1_ * omega(J2_parent_minus_J1_)
              - cross2_ * omega(Jp())
              - trans(J2_parent_minus_J1_));
        } else {
          if (ComputeOrientation) {
            tmpJac_.noalias() = R2().transpose() *
              ( (- cross1_ * R1()) * omega(J1()) + R1() * trans(J1()));
          } else {
            tmpJac_.noalias() =
              ( (- cross1_ * R1()) * omega(J1()) + R1() * trans(J1()));
          }
        }
        jacobian.topRows<3>() = inVel_.rview(tmpJac_);
      }
      if (ComputeOrientation) {
        if (parentJoint_) {
          tmpJac_.noalias() = ( R2().transpose() * Rp() ) * omega(Jp())
            - (R2().transpose() * R1()) * omega(J1());
        } else {
          tmpJac_.noalias() = ( R2().transpose() * R1() ) * omega(J1());
        }
        jacobian.bottomRows<3>() = inVel_.rview(tmpJac_);
      }
    }

    /// Force instanciation of relevant classes
    // template class ExplicitGenericTransformation<               PositionBit | OrientationBit >;
    // template class ExplicitGenericTransformation<               PositionBit                  >;
    // template class ExplicitGenericTransformation<                             OrientationBit >;
    template class ExplicitGenericTransformation< RelativeBit | PositionBit | OrientationBit >;
    template class ExplicitGenericTransformation< RelativeBit | PositionBit                  >;
    // template class ExplicitGenericTransformation< RelativeBit |               OrientationBit >;
  } // namespace constraints
} // namespace hpp
