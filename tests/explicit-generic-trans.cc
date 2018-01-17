// Copyright (c) 2017, Joseph Mirabel
// Authors: Joseph Mirabel (joseph.mirabel@laas.fr)
//
// This file is part of hpp-core.
// hpp-core is free software: you can redistribute it
// and/or modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation, either version
// 3 of the License, or (at your option) any later version.
//
// hpp-core is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Lesser Public License for more details.  You should have
// received a copy of the GNU Lesser General Public License along with
// hpp-core. If not, see <http://www.gnu.org/licenses/>.

#define BOOST_TEST_MODULE ExplicitRelativeTransformationTest
#include <boost/test/included/unit_test.hpp>

#include <boost/shared_ptr.hpp>

// Force benchmark output
#define HPP_ENABLE_BENCHMARK 1
#include <hpp/util/timer.hh>

#include <hpp/constraints/explicit-generic-transformation.hh>

#include <pinocchio/spatial/skew.hpp>
#include <pinocchio/algorithm/joint-configuration.hpp>

#include <hpp/pinocchio/simple-device.hh>
#include <hpp/pinocchio/urdf/util.hh>
#include <hpp/pinocchio/joint.hh>
#include <hpp/pinocchio/configuration.hh>
#include <hpp/pinocchio/liegroup.hh>

// #define PRINT(x) std::cout << x << std::endl
#define PRINT(x) (void)0

using namespace hpp::pinocchio;
using namespace hpp::constraints;

DevicePtr_t createRobot ()
{
  //DevicePtr_t robot = unittest::makeDevice(unittest::HumanoidRomeo, "romeo");
  DevicePtr_t robot (Device::create ("2-objects"));
  urdf::loadModel (robot, 0, "obj1/", "freeflyer", "file://" DATA_DIR "/empty.urdf", "");
  robot->controlComputation((Device::Computation_t) (Device::JOINT_POSITION | Device::JACOBIAN));
  robot->rootJoint()->lowerBound (0, -10);
  robot->rootJoint()->lowerBound (1, -10);
  robot->rootJoint()->lowerBound (2, -10);
  robot->rootJoint()->upperBound (0,  10);
  robot->rootJoint()->upperBound (1,  10);
  robot->rootJoint()->upperBound (2,  10);

  /// Add a freeflyer at the end.
  urdf::loadModel (robot, 0, "obj2/", "freeflyer", "file://" DATA_DIR "/empty.urdf", "");
  JointPtr_t rj = robot->getJointByName("obj2/root_joint");
  rj->lowerBound (0, -10);
  rj->lowerBound (1, -10);
  rj->lowerBound (2, -10);
  rj->upperBound (0,  10);
  rj->upperBound (1,  10);
  rj->upperBound (2,  10);

  return robot;
}

template <typename T> inline typename T::template NRowsBlockXpr<3>::Type trans(const Eigen::MatrixBase<T>& j) { return const_cast<Eigen::MatrixBase<T>&>(j).derived().template topRows   <3>(); }
template <typename T> inline typename T::template NRowsBlockXpr<3>::Type omega(const Eigen::MatrixBase<T>& j) { return const_cast<Eigen::MatrixBase<T>&>(j).derived().template bottomRows<3>(); }

template <typename ExplicitGT> void run_two_freeflyer()
{
  DevicePtr_t robot = createRobot();
  BOOST_REQUIRE (robot);

  JointPtr_t object2 = robot->getJointByName("obj2/root_joint");
  JointPtr_t object1 = robot->getJointByName("obj1/root_joint");

  Transform3f M2inO2 (Transform3f::Identity());
  Transform3f M1inO1 (Transform3f::Identity());

  typename ExplicitGT::Ptr_t ert = ExplicitGT::create (
      "explicit_relative_transformation", robot,
      object1, object2, M1inO1, M2inO2);
  // ExplicitNumericalConstraintPtr_t enm = ert->createNumericalConstraint();

  Configuration_t q     = robot->currentConfiguration (),
                  qrand = se3::randomConfiguration(robot->model()),
                  qout = qrand;

  // Check the output value
  Eigen::RowBlockIndices outConf (ert->outputConf());
  Eigen::RowBlockIndices  inConf (ert-> inputConf());
  Eigen::ColBlockIndices  inVel  (ert-> inputVel ());

  // Compute position of object2 by solving explicit constraints
  LiegroupElement q_obj2 (ert->outputSpace ());
  vector_t q_obj1 (inConf.rview(qout).eval());
  ert->value (q_obj2, q_obj1);
  outConf.lview(qout) = q_obj2.vector ();

  // Test that at solution configuration, object2 and robot frames coincide.
  robot->currentConfiguration(qout);
  robot->computeForwardKinematics();
  Transform3f diff =
    M1inO1.inverse() * object1->currentTransformation().inverse()
    * object2->currentTransformation() * M2inO2;

  if (ExplicitGT::ComputePosition)
    BOOST_CHECK (diff.translation().isZero());
  if (ExplicitGT::ComputeOrientation)
    BOOST_CHECK (diff.rotation().isIdentity());

  // Check Jacobian of implicit numerical constraints by finite difference
  //
  value_type dt (1e-5);
  Configuration_t q0 (qout);
  vector_t v (robot->numberDof ());
  matrix_t J (ert->outputDerivativeSize(), ert->inputDerivativeSize());
  LiegroupElement value0 (ert->outputSpace ()),
    value (ert->outputSpace ());
  ert->value (value0, inConf.rview(q0).eval());
  ert->jacobian (J, inConf.rview(q0).eval());
  PRINT("J=" << std::endl << J);
  PRINT("q0=" << q0.transpose ());
  // First at solution configuration
  for (size_type i=0; i<12; ++i) {
    // test canonical basis vectors
    v.setZero (); v [i] = 1;
    integrate (robot, q0, dt * v, q);
    ert->value (value, inConf.rview(q).eval());
    vector_t df ((value - value0)/dt);
    vector_t Jdq (J * inVel.rviewTranspose(v));
    PRINT("v=" << v.transpose ());
    PRINT("q=" << q.transpose ());
    PRINT("df=" << df.transpose ());
    PRINT("Jdq=" << Jdq.transpose ());
    PRINT("||Jdq - df ||=" << (df - Jdq).norm () << std::endl);
    BOOST_CHECK_SMALL ((df - Jdq).norm (), 1e-4);
    // TODO for ExplicitRelativePosition, the jacobian should be multiplied by R2
    // BOOST_CHECK_SMALL ((object2->currentTransformation().rotation().transpose() * df - Jdq).norm (), 1e-4);
  }
  // Second at random configurations
  for (size_type i=0; i<100; ++i) {
    q0 = se3::randomConfiguration(robot->model());
    ert->jacobian (J, inConf.rview(q0).eval());
    PRINT("J=" << std::endl << J);
    PRINT("q0=" << q0.transpose ());
    ert->value (value0, inConf.rview(q0).eval());
    ert->jacobian (J, inConf.rview(q0).eval());
    for (size_type j=0; j<12; ++j) {
      v.setZero (); v [j] = 1;
      PRINT("v=" << v.transpose ());
      integrate (robot, q0, dt * v, q);
      PRINT("q=" << q.transpose ());
      ert->value (value, inConf.rview(q).eval());
      vector_t df ((value - value0)/dt);
      vector_t Jdq (J * inVel.rviewTranspose(v));
      PRINT("df=" << df.transpose ());
      PRINT("Jdq=" << Jdq.transpose ());
      PRINT("||Jdq - df ||=" << (df - Jdq).norm () << std::endl);
      BOOST_CHECK_SMALL ((df - Jdq).norm (), 1e-4);
    }
    PRINT("");
  }
}

BOOST_AUTO_TEST_CASE (two_freeflyer_trans) { run_two_freeflyer<ExplicitRelativeTransformation>(); }
BOOST_AUTO_TEST_CASE (two_freeflyer_pos  ) { run_two_freeflyer<ExplicitRelativePosition      >(); }

template <typename ExplicitGT> void run_two_frames_on_freeflyer()
{
  DevicePtr_t robot = createRobot();
  BOOST_REQUIRE (robot);

  JointPtr_t object2 = robot->getJointByName("obj2/root_joint");
  JointPtr_t object1  = robot->getJointByName("obj1/root_joint");

  Transform3f M2inO2 (Transform3f::Random ());
  Transform3f M1inO1 (Transform3f::Random ());

  PRINT("M2inO2=" << M2inO2);
  PRINT("M1inO1=" << M1inO1);

  typename ExplicitGT::Ptr_t ert = ExplicitGT::create (
      "explicit_relative_transformation", robot,
      object1, object2, M1inO1, M2inO2);
  // ExplicitNumericalConstraintPtr_t enm = ert->createNumericalConstraint();

  Configuration_t q     = robot->currentConfiguration (),
                  qrand = se3::randomConfiguration(robot->model()),
                  qout = qrand;

  // Check the output value
  Eigen::RowBlockIndices outConf (ert->outputConf());
  Eigen::RowBlockIndices  inConf (ert-> inputConf());
  Eigen::ColBlockIndices  inVel  (ert-> inputVel ());

  // Compute position of object2 by solving explicit constraints
  LiegroupElement q_obj2 (ert->outputSpace ());
  vector_t q_obj1 (inConf.rview(qout).eval());
  ert->value (q_obj2, q_obj1);
  outConf.lview(qout) = q_obj2.vector ();

  // Test that at solution configuration, object2 and robot frames coincide.
  robot->currentConfiguration(qout);
  robot->computeForwardKinematics();
  Transform3f diff =
    M1inO1.inverse() * object1->currentTransformation().inverse()
    * object2->currentTransformation() * M2inO2;

  if (ExplicitGT::ComputePosition) {
    BOOST_CHECK (diff.translation().isZero());
    vector3_t diff2 = object1->currentTransformation().rotation() * M1inO1.translation()
      - (object2->currentTransformation().translation() + M2inO2.translation());
    BOOST_CHECK (diff2.isZero());
  }
  if (ExplicitGT::ComputeOrientation)
    BOOST_CHECK (diff.rotation().isIdentity());

  // Check Jacobian of implicit numerical constraints by finite difference
  //
  value_type dt (1e-5);
  Configuration_t q0 (qout);
  vector_t v (robot->numberDof ());
  matrix_t J (ert->outputDerivativeSize(), ert->inputDerivativeSize());
  LiegroupElement value0 (ert->outputSpace ()),
    value (ert->outputSpace ());
  ert->value (value0, inConf.rview(q0).eval());
  ert->jacobian (J, inConf.rview(q0).eval());
    PRINT("J=" << std::endl << J);
    PRINT("q0=" << q0.transpose ());
  // First at solution configuration
  for (size_type i=0; i<12; ++i) {
    // test canonical basis vectors
    v.setZero (); v [i] = 1;
    integrate (robot, q0, dt * v, q);
    ert->value (value, inConf.rview(q).eval());
    vector_t df ((value - value0)/dt);
    vector_t Jdq (J * inVel.rviewTranspose(v));
    PRINT("v=" << v.transpose ());
    PRINT("q=" << q.transpose ());
    PRINT("df=" << df.transpose ());
    PRINT("Jdq=" << Jdq.transpose ());
    PRINT("||Jdq - df ||=" << (df - Jdq).norm () << std::endl);
    BOOST_CHECK_SMALL ((df - Jdq).norm (), 1e-4);
  }
  // Second at random configurations
  for (size_type i=0; i<100; ++i) {
    q0 = se3::randomConfiguration(robot->model());
    ert->value (value0, inConf.rview(q0).eval());
    ert->jacobian (J, inConf.rview(q0).eval());
    PRINT("J=" << std::endl << J);
    PRINT("q0=" << q0.transpose ());
    for (size_type j=0; j<12; ++j) {
      v.setZero (); v [j] = 1;
      PRINT("v=" << v.transpose ());
      integrate (robot, q0, dt * v, q);
      PRINT("q=" << q.transpose ());
      ert->value (value, inConf.rview(q).eval());
      vector_t df ((value - value0)/dt);
      vector_t Jdq (J * inVel.rviewTranspose(v));
      PRINT("df=" << df.transpose ());
      PRINT("Jdq=" << Jdq.transpose ());
      PRINT("||Jdq - df ||=" << (df - Jdq).norm () << std::endl);
      BOOST_CHECK_SMALL ((df - Jdq).norm (), 1e-4);
    }
    PRINT("");
  }
}

BOOST_AUTO_TEST_CASE (two_frames_on_freeflyer_trans) { run_two_frames_on_freeflyer<ExplicitRelativeTransformation>(); }
BOOST_AUTO_TEST_CASE (two_frames_on_freeflyer_pos  ) { run_two_frames_on_freeflyer<ExplicitRelativePosition      >(); }
