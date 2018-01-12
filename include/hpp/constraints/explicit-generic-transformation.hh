//
// Copyright (c) 2016 CNRS
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

#ifndef HPP_CONSTRAINTS_EXPLICIT_GENERIC_TRANSFORMATION_HH
# define HPP_CONSTRAINTS_EXPLICIT_GENERIC_TRANSFORMATION_HH

# include <pinocchio/spatial/se3.hpp>

# include <hpp/pinocchio/joint.hh>

# include <hpp/constraints/fwd.hh>
# include <hpp/constraints/config.hh>
# include <hpp/constraints/differentiable-function.hh>
# include <hpp/constraints/matrix-view.hh>

namespace hpp {
  namespace constraints {
    /// \addtogroup constraints
    /// \{

    /** ExplicitGenericTransformation class encapsulates 6 possible differentiable
     *  functions: ExplicitPosition, ExplicitOrientation, ExplicitTransformation
     *  and their relative versions ExplicitRelativePosition,
     *  ExplicitRelativeOrientation, ExplicitRelativeTransformation.
     *  At the time of writting, only the relative functions are supported.
     *
     *  These functions compute the position of frame \ref frame2InJoint2
     *  in joint \ref joint2 frame, in the frame \ref frame1InJoint1
     *  in \ref joint1 frame. For absolute functions, \ref joint1 is
     *  NULL and joint1 frame is the world frame.
     *
     *  The value of the ExplicitRelativeTransformation function is an element
     *  of SE(3) corresponding to the configuration parameter of \ref joint2
     *
    */
    template <int _Options>
    class HPP_CONSTRAINTS_DLLAPI ExplicitGenericTransformation :
      public DifferentiableFunction
    {
    public:
      typedef boost::shared_ptr <ExplicitGenericTransformation> Ptr_t;
      typedef boost::weak_ptr <ExplicitGenericTransformation> WkPtr_t;

      enum {
        // IsRelative         = _Options & RelativeBit,
        IsRelative         = 1,
        ComputeOrientation = _Options & OrientationBit,
        ComputePosition    = _Options & PositionBit,
        IsPosition         = ComputePosition  && !ComputeOrientation,
        IsOrientation      = !ComputePosition && ComputeOrientation,
        IsTransform        = ComputePosition  && ComputeOrientation,
        MaskSize           = (ComputePosition?3:0)
      };

      /// \cond
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW
      /// \endcond

      /// Object builder for relative functions.
      ///
      /// \param name the name of the constraints,
      /// \param robot the robot the constraints is applied to,
      /// \param joint1 the first joint the transformation of which is
      ///               constrained,
      /// \param joint2 the second joint the transformation of which is
      ///               constrained,
      /// \param reference desired relative transformation
      ///        \f$T_1(\mathbf{q})^{-1} T_2(\mathbf{q})\f$ between the joints.
      /// \param mask which component of the error vector to take into
      ///        account.
      static Ptr_t create (const std::string& name, const DevicePtr_t& robot,
          const JointConstPtr_t& joint1, const JointConstPtr_t& joint2,
          const Transform3f& reference,
          std::vector <bool> mask = std::vector<bool>(MaskSize,true));

      /// Object builder for relative functions.
      ///
      /// \param name the name of the constraints,
      /// \param robot the robot the constraints is applied to,
      /// \param joint1 the first joint the transformation of which is
      ///               constrained,
      /// \param joint2 the second joint the transformation of which is
      ///               constrained,
      /// \param frame1 position of a fixed frame in joint 1,
      /// \param frame2 position of a fixed frame in joint 2,
      /// \param mask vector of 6 boolean defining which coordinates of the
      ///        error vector to take into account.
      /// \note if joint1 is 0x0, joint 1 frame is considered to be the global
      ///       frame.
      ///
      /// \note For RelativePosition, the rotation part of frame1 defines the
      ///       frame in which the error is expressed and the rotation of frame2
      ///       has no effect.
      static Ptr_t create (const std::string& name, const DevicePtr_t& robot,
	 const JointConstPtr_t& joint1,  const JointConstPtr_t& joint2,
	 const Transform3f& frame1, const Transform3f& frame2,
         std::vector <bool> mask = std::vector<bool>(MaskSize,true));

      virtual ~ExplicitGenericTransformation () {}

      /// Set desired relative transformation of joint2 in joint1
      ///
      inline void reference (const Transform3f& reference)
      {
        ref_ = reference;
      }

      /// Get desired relative orientation
      inline const Transform3f& reference () const
      {
	return ref_;
      }

      const segments_t& inputConf () const
      {
        return inConf_.indices();
      }
      const segments_t& inputVel () const
      {
        return inVel_.indices();
      }
      const segments_t& outputConf () const
      {
        return outConf_.indices();
      }
      const segments_t& outputVel () const
      {
        return outVel_.indices();
      }

      virtual std::ostream& print (std::ostream& o) const;

    protected:
      ///Constructor
      ///
      /// \param name the name of the constraints,
      /// \param robot the robot the constraints is applied to,
      ///        \f$T_1(\mathbf{q})^{-1} T_2(\mathbf{q})\f$ between the joints.
      /// \param mask vector of 6 boolean defining which coordinates of the
      ///        error vector to take into account.
      ExplicitGenericTransformation (
          const std::string& name      , const DevicePtr_t& robot,
          const JointConstPtr_t& joint1, const JointConstPtr_t& joint2,
          const Transform3f& reference,
          const segments_t inConf , const segments_t outConf,
          const segments_t inVel  , const segments_t outVel ,
          std::vector <bool> mask);

      void init (const WkPtr_t& self)
      {
        self_ = self;
      }

      /// Compute value of error
      ///
      /// \param argument configuration of the robot,
      /// \retval result error vector
      virtual void impl_compute	(LiegroupElement& result,
				 ConfigurationIn_t argument) const throw ();
      virtual void impl_jacobian (matrixOut_t jacobian,
				  ConfigurationIn_t arg) const throw ();
    private:
      void forwardKinematics (const vectorIn_t& arg) const;
      void computeError (const ConfigurationIn_t& argument) const;

      inline const Transform3f& M1 () const { return joint1_->currentTransformation(); }
      inline const matrix3_t& R1 () const { return M1().rotation(); }
      inline const vector3_t& t1 () const { return M1().translation(); }
      inline const JointJacobian_t& J1 () const { return joint1_->jacobian(); }

      inline const Transform3f& M2 () const { return joint2_->currentTransformation(); }
      inline const vector3_t& t2 () const { return M2().translation(); }
      inline const matrix3_t& R2 () const { return M2().rotation(); }
      inline const JointJacobian_t& J2 () const { return joint2_->jacobian(); }

      inline const Transform3f& Mp () const { return parentJoint_->currentTransformation(); }
      inline const vector3_t& tp () const { return Mp().translation(); }
      inline const matrix3_t& Rp () const { return Mp().rotation(); }
      inline const JointJacobian_t& Jp () const { return parentJoint_->jacobian(); }

      inline const Transform3f& pM2 () const { return joint2_->positionInParentFrame(); }
      inline const vector3_t& pt2 () const { return pM2().translation(); }
      inline const matrix3_t& pR2 () const { return pM2().rotation(); }

      DevicePtr_t robot_;
      JointConstPtr_t joint1_, joint2_, parentJoint_;
      Transform3f ref_; // F1inJ1 * F2inJ2.inverse()

      Eigen::RowBlockIndices inConf_;
      Eigen::ColBlockIndices inVel_;
      Eigen::RowBlockIndices outConf_, outVel_;

      const std::vector<bool> mask_;
      WkPtr_t self_;
      mutable vector_t qsmall_, q_;
      mutable matrix3_t cross1_, cross2_;
      mutable JointJacobian_t J2_parent_minus_J1_;
      mutable Eigen::Matrix<value_type, 3, Eigen::Dynamic> tmpJac_;
      mutable Transform3f freeflyerPose_;
    }; // class GenericTransformation
    /// \}
  } // namespace constraints
} // namespace hpp
#endif // HPP_CONSTRAINTS_EXPLICIT_GENERIC_TRANSFORMATION_HH
