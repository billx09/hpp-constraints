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

#include <hpp/constraints/generic-transformation.hh>

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
      /// ------- Generic Transform Data ---------------------------------------
      template <bool rel> struct GTOriDataV {};
      template <> struct GTOriDataV<true>
      {
        Transform3f M;
        value_type theta;
      };
      template <bool rel> struct GTOriDataJ {};
      template <> struct GTOriDataJ<true>
      {
        eigen::matrix3_t JlogXTR1inJ1;
      };
      /// This class contains the data of the GenericTransformation class.
      template <bool rel> struct GTDataBase
      {
        hpp::pinocchio::DeviceSync device;
        const GenericTransformationModel<rel>& model;
        hpp::pinocchio::DeviceData& ddata () { return device.d(); }

        const JointJacobian_t& J2 () { return model.joint2->jacobian(ddata()); }
        const Transform3f& M2 () { return model.joint2->currentTransformation(ddata()); }
        const vector3_t& t2 () { return model.joint2->currentTransformation(ddata()).translation(); }
        const matrix3_t& R2 () { return model.joint2->currentTransformation(ddata()).rotation(); }

        const JointJacobian_t& J1 () { return model.getJoint1()->jacobian(this->ddata()); }
        const Transform3f& M1 () { return model.getJoint1()->currentTransformation(this->ddata()); }
        const matrix3_t& R1 () { return model.getJoint1()->currentTransformation(this->ddata()).rotation(); }
        const vector3_t& t1 () { return model.getJoint1()->currentTransformation(this->ddata()).translation(); }

        GTDataBase (const GenericTransformationModel<rel>& m, const DevicePtr_t& d)
          : device (d), model(m) {}
      };
      template <bool rel, bool pos, bool ori> struct GTDataV :
        GTDataBase<rel>, GTOriDataV<ori>
      {
        enum {
          NbRows = (pos?3:0)+(ori?3:0)
        };
        typedef Eigen::Matrix<value_type, NbRows, 1> ValueType;
        ValueType value;

        GTDataV (const GenericTransformationModel<rel>& m, const DevicePtr_t& d)
          : GTDataBase<rel>(m, d) {}
      };
      /// This class contains the data of the GenericTransformation class.
      template <bool rel, bool pos, bool ori> struct GTDataJ :
        GTDataV<rel, pos, ori>,
        GTOriDataJ<ori>
      {
        enum {
          NbRows = (pos?3:0)+(ori?3:0),
          RowPos = (pos? 0:-1),
          RowOri = (ori?(pos?3:0):-1)
        };
        typedef Eigen::Matrix<value_type, NbRows, Eigen::Dynamic> JacobianType;
        JacobianType jacobian;
        Eigen::Matrix<value_type, 3, Eigen::Dynamic> tmpJac;
        eigen::vector3_t cross1, cross2;

        GTDataJ (const GenericTransformationModel<rel>& m, const DevicePtr_t& d)
          : GTDataV<rel,pos,ori> (m, d)
        {
          if (!m.fullPos || !m.fullOri) jacobian.resize ((int)NbRows, m.cols);
          cross1.setZero();
          if (m.t2isZero) cross2.setZero();
        }
      };

      /// ------- Compute log --------------------------------------------------
      /** Compute jacobian of function log of rotation matrix in SO(3)

          Let us consider a matrix
          \f$R=\exp \left[\mathbf{r}\right]_{\times}\in SO(3)\f$.
          This functions computes the Jacobian of the function from
          \f$SO(3)\f$ into \f$\mathbf{R}^3\f$ that maps \f$R\f$ to
          \f$\mathbf{r}\f$. In other words,
          \f{equation*}
          \dot{\mathbf{r}} = J_{log}(R)\ \omega\,\,\,\mbox{with}\,\,\,
          \dot {R} = \left[\omega\right]_{\times} R
          \f}
          \warning Two representations of the angular velocity \f$\omega\f$ are
                   possible:
                   \li \f$\dot{R} = \left[\omega\right]_{\times}R\f$ or
                   \li \f$\dot{R} = R\left[\omega\right]_{\times}\f$.

                   The expression below assumes the first representation is
                   used.
          \param theta angle of rotation \f$R\f$, also \f$\|r\|\f$,
          \param log 3d vector \f$\mathbf{r}\f$,
          \retval Jlog matrix \f$J_{log} (R)\f$.

          \f{align*}
          J_{log} (R) &=& \frac{\|\mathbf{r}\|\sin\|\mathbf{r}\|}{2(1-\cos\|\mathbf{r}\|)} I_3 - \frac {1}{2}\left[\mathbf{r}\right]_{\times} + (\frac{1}{\|\mathbf{r}\|^2} - \frac{\sin\|\mathbf{r}\|}{2\|\mathbf{r}\|(1-\cos\|\mathbf{r}\|)}) \mathbf{r}\mathbf{r}^T\\
           &=& I_3 -\frac{1}{2}\left[\mathbf{r}\right]_{\times} +  \left(\frac{2(1-\cos\|\mathbf{r}\|) - \|\mathbf{r}\|\sin\|\mathbf{r}\|}{2\|\mathbf{r}\|^2(1-\cos\|\mathbf{r}\|)}\right)\left[\mathbf{r}\right]_{\times}^2
           \f} */
      template <typename Derived>
      void computeJlog (const value_type& theta, const Eigen::MatrixBase<Derived>& log, matrix3_t& Jlog)
      {
        if (theta < 1e-6)
          Jlog.setIdentity();
        else {
          // Jlog = alpha I
          const value_type ct = cos(theta), st = sin(theta);
          const value_type st_1mct = st/(1-ct);

          Jlog.setZero ();
          Jlog.diagonal().setConstant (theta*st_1mct);

          // Jlog += -r_{\times}/2
          Jlog(0,1) =  log(2); Jlog(1,0) = -log(2);
          Jlog(0,2) = -log(1); Jlog(2,0) =  log(1);
          Jlog(1,2) =  log(0); Jlog(2,1) = -log(0);
          Jlog /= 2;

          const value_type alpha = 1/(theta*theta) - st_1mct/(2*theta);
          Jlog.noalias() += alpha * log * log.transpose ();
        }
      }

      typedef JointJacobian_t::ConstNRowsBlockXpr<3>::Type HalfJacobian_t;
      inline HalfJacobian_t omega(const JointJacobian_t& j) { return j.bottomRows<3>(); }
      inline HalfJacobian_t trans(const JointJacobian_t& j) { return j.topRows<3>(); }

      static inline size_type size (std::vector<bool> mask)
      {
        size_type res = 0;
        for (std::vector<bool>::iterator it = mask.begin (); it != mask.end ();
            ++it) {
          if (*it) ++res;
        }
        return res;
      }

      template <bool flag /* false */ > struct unary
      {
        template <bool rel, bool pos> static inline void log (
            GTDataV<rel, pos, flag>&) {}
        template <bool rel, bool pos> static inline void Jlog (
            GTDataJ<rel, pos, flag>&) {}
      };
      template <> struct unary <true>
      {
        template <bool rel, bool pos> static inline void log (
            GTDataV<rel, pos, true>& d)
          {
            logSO3(d.M.rotation(), d.theta, d.value.template tail<3>());
            hppDnum (info, "theta=" << d.theta);
          }
        template <bool rel, bool pos> static inline void Jlog (
            GTDataJ<rel, pos, true>& d)
          {
            computeJlog(d.theta, d.value.template tail<3>(), d.JlogXTR1inJ1);
            hppDnum (info, "Jlog_: " << d.JlogXTR1inJ1);
            if (!d.model.R1isID) d.JlogXTR1inJ1 *= d.model.F1inJ1.rotation().transpose();
          }
      };

      template <bool ori, typename Data, typename Derived> void assign_if
        (bool cond, Data& d, matrixOut_t J,
         const Eigen::MatrixBase<Derived>& rhs,
         const size_type& startRow)
      {
        const int& rowCache = (ori ? Data::RowOri : Data::RowPos);
        if (cond) d.jacobian.template middleRows<3>(rowCache)                       .noalias() = rhs;
        else               J.template middleRows<3>(startRow).leftCols(d.model.cols).noalias() = rhs;
      }

      /// ------- Compute jacobian ---------------------------------------------
      template <bool lflag /*rel*/, bool rflag /*false*/> struct binary
      {
        // the first template allow us to consider relative transformation as
        // absolute when joint1 is NULL, at run time
        template <bool rel, bool pos> static inline void Jorientation (
            GTDataJ<rel, pos, rflag>&, matrixOut_t) {}
        template <bool rel, bool ori> static inline void Jtranslation (
            GTDataJ<rel, rflag, ori>&, matrixOut_t) {}
      };
      template <> struct binary<false, true> // Absolute
      {
        template <bool rel, bool pos> static inline void Jorientation (
            GTDataJ<rel, pos, true>& d, matrixOut_t J)
        {
          assign_if<true>(!d.model.fullOri, d, J,
            (d.JlogXTR1inJ1 * d.R2()) * omega(d.J2()),
            d.model.rowOri);
        }
        template <bool rel, bool ori> static inline void Jtranslation (
            GTDataJ<rel, true, ori>& d,
            matrixOut_t J)
        {
          const JointJacobian_t& J2 (d.J2());
          const matrix3_t& R2 (d.R2());
          const matrix3_t& R1inJ1 (d.model.F1inJ1.rotation ());

          // hpp-model: J = 1RT* ( 0Jt2 - [ 0R2 2t* ]x 0Jw2 )
          // pinocchio: J = 1RT* ( 0R2 2Jt2 - [ 0R2 2t* ]x 0R2 2Jw2 )
          if (!d.model.t2isZero) {
            d.tmpJac.noalias() = ( R2.colwise().cross(d.cross2)) * omega(J2);
            d.tmpJac.noalias() += R2 * trans(J2);
            if (d.model.R1isID) {
              assign_if<false> (!d.model.fullPos, d, J, d.tmpJac, 0);
            } else { // Generic case
              assign_if<false> (!d.model.fullPos, d, J, R1inJ1.transpose() * d.tmpJac, 0);
            }
          } else {
            if (d.model.R1isID)
              assign_if<false> (!d.model.fullPos, d, J, R2 * trans(J2), 0);
            else
              assign_if<false> (!d.model.fullPos, d, J, (R1inJ1.transpose() * R2) * trans(J2), 0);
          }
        }
      };
      template <> struct binary<true, true> // Relative
      {
        template <bool pos> static inline void Jorientation (
            GTDataJ<true, pos, true>& d,
            matrixOut_t J)
        {
          d.tmpJac.noalias() = (d.R1().transpose() * d.R2()) * omega(d.J2());
          d.tmpJac.noalias() -= omega(d.J1());
          assign_if<true>(!d.model.fullOri, d, J,
              d.JlogXTR1inJ1 * d.tmpJac,
              d.model.rowOri);
        }
        template <bool ori> static inline void Jtranslation (
            GTDataJ<true, true, ori>& d,
            matrixOut_t J)
        {
          const JointJacobian_t& J1 (d.J1()); const JointJacobian_t& J2 (d.J2());
          const matrix3_t&       R1 (d.R1()); const matrix3_t&       R2 (d.R2());
          const matrix3_t& R1inJ1 (d.model.F1inJ1.rotation ());

          // J = 1RT* 0RT1 ( A + B )
          // hpp-model:
          // A = [ 0t2 - 0t1 0R2 2t* ]x 0Jw1
          // B = ( 0Jt2 - 0Jt1 - [ 0R2 2t* ]x 0Jw2 )
          // pinocchio:
          // A = [ 0t2 - 0t1 0R2 2t* ]x 0R1 1Jw1
          // B = ( 0R2 2Jt2 - 0R1 1Jt1 - [ 0R2 2t* ]x 0R2 2Jw2 )
          d.tmpJac.noalias() = (- R1.transpose() * R1.colwise().cross(d.cross1) ) * omega(J1); // A
          d.tmpJac.noalias() += ( R1.transpose() * R2 ) * trans(J2);  // B1
          d.tmpJac.noalias() -= trans(J1); // B2
          if (!d.model.t2isZero)
            d.tmpJac.noalias() += R1.transpose() * R2.colwise().cross(d.cross2) * omega(J2); // B3
          if (d.model.R1isID) assign_if<false>(!d.model.fullPos, d, J,                      d.tmpJac, 0);
          else                assign_if<false>(!d.model.fullPos, d, J, R1inJ1.transpose() * d.tmpJac, 0);
        }
      };

      /// ------- Compute relative transform -----------------------------------
      template <bool compileTimeRel /* false */, bool ori /* false */> struct relativeTransform {
        template <bool runtimeRel> static inline void run (GTDataV<runtimeRel, true, false>& d)
        {
          // There is no joint1
          const Transform3f& M2 = d.M2 ();
          d.value.noalias() = M2.act (d.model.F2inJ2.translation());
          if (!d.model.t1isZero) d.value.noalias() -= d.model.F1inJ1.translation();
          if (!d.model.R1isID)
            d.value.applyOnTheLeft(d.model.F1inJ1.rotation().transpose());
        }
      };
      template <> struct relativeTransform<false, true> {
        template <bool runtimeRel, bool pos> static inline void run (
            GTDataV<runtimeRel, pos, true>& d)
        {
          const Transform3f& M2 = d.M2 ();
          d.M = d.model.F1inJ1.actInv(M2 * d.model.F2inJ2);
          if (pos) d.value.template head<3>().noalias() = d.M.translation();
        }
      };
      template <> struct relativeTransform<true, true> {
        template <bool pos> static inline void run (
            GTDataV<true, pos, true>& d)
        {
          if (d.model.joint1 == NULL) {
            // runtime absolute reference.
            relativeTransform<false, true>::run(d);
            return;
          }
          const Transform3f& M1 = d.M1 ();
          const Transform3f& M2 = d.M2 ();
          d.M = d.model.F1inJ1.actInv(M1.actInv(M2 * d.model.F2inJ2));
          if (pos) d.value.template head<3>().noalias() = d.M.translation();
        }
      };
      template <> struct relativeTransform<true, false> {
        static inline void run (GTDataV<true, true, false>& d)
        {
          if (d.model.joint1 == NULL) {
            // runtime absolute reference.
            relativeTransform<false, false>::run(d);
            return;
          }
          const Transform3f& M2 = d.M2 ();
          const Transform3f& M1 = d.M1 ();
          d.value.noalias() = M2.act (d.model.F2inJ2.translation())
                              - M1.translation();
          d.value.applyOnTheLeft(M1.rotation().transpose());

          if (!d.model.t1isZero) d.value.noalias() -= d.model.F1inJ1.translation();
          if (!d.model.R1isID)
            d.value.applyOnTheLeft(d.model.F1inJ1.rotation().transpose());
        }
      };

      template <bool rel, bool pos, bool ori> struct compute
      {
        static inline void error (GTDataV<rel, pos, ori>& d)
        {
          relativeTransform<rel, ori>::run (d);
          unary<ori>::log(d);
        }

        static inline void jacobian (GTDataJ<rel, pos, ori>& d,
            matrixOut_t jacobian, const std::vector<bool>& mask)
        {
          const Transform3f& M2 = d.M2 ();
          const vector3_t& t2inJ2 (d.model.F2inJ2.translation ());
          const vector3_t& t2 (M2.translation ());
          const matrix3_t& R2 (M2.rotation ());

          if (!d.model.t2isZero)
            d.cross2.noalias() = R2*t2inJ2;

          unary<ori>::Jlog (d);

          // rel:           relative known at compile time
          // d.getJoint1(): relative known at run time
          if (rel && d.model.getJoint1()) {
            const Transform3f& M1 = d.M1();
            const vector3_t& t1 (M1.translation ());
            d.cross1.noalias() = d.cross2 + t2 - t1;
            binary<rel, pos>::Jtranslation (d, jacobian);
            binary<rel, ori>::Jorientation (d, jacobian);
          } else {
            d.cross1.noalias() = d.cross2 + t2;
            binary<false, pos>::Jtranslation (d, jacobian);
            binary<false, ori>::Jorientation (d, jacobian);
          }

          // Copy necessary rows.
          size_type index=0;
          const size_type lPos = (pos?3:0), lOri = (ori?3:0);
          if (!d.model.fullPos) {
            for (size_type i=0; i<lPos; ++i) {
              if (mask [i]) {
                jacobian.row(index).leftCols(d.model.cols).noalias() = d.jacobian.row(i); ++index;
              }
            }
          } else index = lPos;
          if (!d.model.fullOri) {
            for (size_type i=lPos; i<lPos+lOri; ++i) {
              if (mask [i]) {
                jacobian.row(index).leftCols(d.model.cols).noalias() = d.jacobian.row(i); ++index;
              }
            }
          }
          jacobian.rightCols(jacobian.cols()-d.model.cols).setZero();
        }
      };
    }

    template <int _Options> std::ostream&
      GenericTransformation<_Options>::print (std::ostream& os) const
    {
      os << (IsRelative ? "Relative" : "") <<
            (ComputePosition ? (ComputeOrientation ? "Transformation" : "Position")
                   : "Orientation") << ": " << name()
        << ", active dof: "
        << pretty_print (BlockIndex::fromLogicalExpression (activeParameters_)) << incindent
        << iendl << "Joint1: "         << ((IsRelative && joint1()) ? joint1()->name() : "World")
        << iendl << "Frame in joint 1" << incindent << iendl; pinocchio::display(os, frame1InJoint1()) << decindent
        << iendl << "Joint2: "         << joint2()->name()
        << iendl << "Frame in joint 2" << incindent << iendl; pinocchio::display(os, frame2InJoint2()) << decindent
        << iendl << "mask: ";
      for (size_type i=0; i<ValueSize; ++i) os << mask_ [i] << ", ";
      return os << decindent;
    }

    template <int _Options> typename GenericTransformation<_Options>::Ptr_t
      GenericTransformation<_Options>::create
    (const std::string& name, const DevicePtr_t& robot,
     const JointConstPtr_t& joint2,
     const Transform3f& reference, std::vector <bool> mask)
    {
      GenericTransformation<_Options>* ptr =
        new GenericTransformation<_Options> (name, robot, mask);
      ptr->joint1 (JointConstPtr_t());
      ptr->joint2 (joint2);
      ptr->reference (reference);
      Ptr_t shPtr (ptr);
      ptr->init (shPtr);
      return shPtr;
    }

    template <int _Options> typename GenericTransformation<_Options>::Ptr_t
      GenericTransformation<_Options>::create
    (const std::string& name, const DevicePtr_t& robot,
     /* World frame          */ const JointConstPtr_t& joint2,
     const Transform3f& frame2, const Transform3f& frame1,
     std::vector <bool> mask)
    {
      GenericTransformation<_Options>* ptr =
        new GenericTransformation<_Options> (name, robot, mask);
      ptr->joint1 (JointConstPtr_t());
      ptr->joint2 (joint2);
      ptr->frame1InJoint1 (frame1);
      ptr->frame2InJoint2 (frame2);
      Ptr_t shPtr (ptr);
      return shPtr;
    }

    template <int _Options> typename GenericTransformation<_Options>::Ptr_t
      GenericTransformation<_Options>::create
    (const std::string& name, const DevicePtr_t& robot,
     const JointConstPtr_t& joint1, const JointConstPtr_t& joint2,
     const Transform3f& reference, std::vector <bool> mask)
    {
      GenericTransformation<_Options>* ptr =
        new GenericTransformation<_Options> (name, robot, mask);
      ptr->joint1 (joint1);
      ptr->joint2 (joint2);
      ptr->reference (reference);
      Ptr_t shPtr (ptr);
      ptr->init (shPtr);
      return shPtr;
    }

    template <int _Options> typename GenericTransformation<_Options>::Ptr_t
      GenericTransformation<_Options>::create
    (const std::string& name, const DevicePtr_t& robot,
     const JointConstPtr_t& joint1, const JointConstPtr_t& joint2,
     const Transform3f& frame1, const Transform3f& frame2,
     std::vector <bool> mask)
    {
      GenericTransformation<_Options>* ptr =
        new GenericTransformation<_Options> (name, robot, mask);
      ptr->joint1 (joint1);
      ptr->joint2 (joint2);
      ptr->frame1InJoint1 (frame1);
      ptr->frame2InJoint2 (frame2);
      Ptr_t shPtr (ptr);
      return shPtr;
    }

    template <int _Options>
      GenericTransformation<_Options>::GenericTransformation
      (const std::string& name, const DevicePtr_t& robot,
       std::vector <bool> mask) :
        DifferentiableFunction (robot->configSize (), robot->numberDof (),
                                LiegroupSpace::Rn (size (mask)), name),
        robot_ (robot),
        m_(robot->numberDof()-robot->extraConfigSpace().dimension()),
        mask_ (mask)
    {
      assert(mask.size()==ValueSize);
      std::size_t iOri = 0;
      m_.rowOri = 0;
      if (ComputePosition) {
        for (size_type i=0; i<3; ++i) if (mask_[i]) m_.rowOri++;
        m_.fullPos = (m_.rowOri==3);
        iOri = 3;
      } else m_.fullPos = false;
      if (ComputeOrientation)
        m_.fullOri = mask_[iOri + 0] && mask_[iOri + 1] && mask_[iOri + 2];
      else m_.fullOri = false;
    }

    template <int _Options>
    inline void GenericTransformation<_Options>::computeActiveParams ()
    {
      typedef se3::JointIndex JointIndex;
      activeParameters_.setConstant (false);
      activeDerivativeParameters_.setConstant (false);

      const se3::Model& model = robot_->model();
      const JointIndex id1 = (joint1() ? joint1()->index() : 0),
                       id2 = (joint2() ? joint2()->index() : 0);
      JointIndex i1 = id1, i2 = id2;

      std::vector<JointIndex> from1, from2;
      while (i1 != i2)
      {
        JointIndex i;
        if (i1 > i2) {
          i = i1;
          i1 = model.parents[i1];
        } else /* if (i1 < i2) */ {
          i = i2;
          i2 = model.parents[i2];
        }
        if (i > 0) {
          activeParameters_
            .segment(model.joints[i].idx_q(), 
                model.joints[i].nq()
                ).setConstant(true);
          activeDerivativeParameters_
            .segment(model.joints[i].idx_v(), 
                model.joints[i].nv()
                ).setConstant(true);
        }
      }
      assert (i1 == i2);
    }

    template <int _Options>
    void GenericTransformation<_Options>::impl_compute
    (LiegroupElement& result, ConfigurationIn_t argument) const throw ()
    {
      GTDataV<IsRelative, ComputePosition, ComputeOrientation> data (m_, robot_);

      data.device.currentConfiguration (argument);
      data.device.computeForwardKinematics ();
      compute<IsRelative, ComputePosition, ComputeOrientation>::error (data);

      size_type index=0;
      for (size_type i=0; i<ValueSize; ++i) {
	if (mask_ [i]) {
	  result.vector () [index] = data.value[i]; ++index;
	}
      }
    }

    template <int _Options>
    void GenericTransformation<_Options>::impl_jacobian
    (matrixOut_t jacobian, ConfigurationIn_t arg) const throw ()
    {
      // TODO there is still a little bit a memory which is dynamically
      // allocated in GTDataJ. At the moment, this allocation is necessary to
      // support multithreadind. To avoid it, DeviceData should provide some
      // a temporary buffer to pass to an Eigen::Map
      GTDataJ<IsRelative, ComputePosition, ComputeOrientation> data (m_, robot_);

      data.device.currentConfiguration (arg);
      data.device.computeForwardKinematics ();
      compute<IsRelative, ComputePosition, ComputeOrientation>::error (data);
      compute<IsRelative, ComputePosition, ComputeOrientation>::jacobian (data, jacobian, mask_);

#ifdef CHECK_JACOBIANS
      const value_type eps = std::sqrt(Eigen::NumTraits<value_type>::epsilon());
      matrix_t Jfd (outputDerivativeSize(), inputDerivativeSize());
      Jfd.setZero();
      finiteDifferenceCentral(Jfd, arg, robot_, eps);
      size_type row, col;
      value_type maxError = (jacobian - Jfd).cwiseAbs().maxCoeff(&row,&col);
      if (maxError > std::sqrt(eps)) {
        hppDout (error, "Jacobian of " << name() << " does not match central finite difference. "
            "DOF " << col << " at row " << row << ": "
            << maxError << " > " << /* HESSIAN_MAXIMUM_COEF << " * " << */ std::sqrt(eps)
            );
        hppDnum (error,
            "Jacobian is" << iendl << jacobian << iendl
            "Finite diff is" << iendl << Jfd << iendl
            "Difference is" << iendl << (jacobian - Jfd));
      }
#endif
    }

    /// Force instanciation of relevant classes
    template class GenericTransformation<               PositionBit | OrientationBit >;
    template class GenericTransformation<               PositionBit                  >;
    template class GenericTransformation<                             OrientationBit >;
    template class GenericTransformation< RelativeBit | PositionBit | OrientationBit >;
    template class GenericTransformation< RelativeBit | PositionBit                  >;
    template class GenericTransformation< RelativeBit |               OrientationBit >;
  } // namespace constraints
} // namespace hpp
