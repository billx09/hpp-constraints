#include <hpp/constraints/convex-collision-avoidance.hh>

#include <iostream>

#include <pinocchio/multibody/model.hpp>
#include <pinocchio/multibody/geometry.hpp>
#include <pinocchio/algorithm/joint-configuration.hpp>

#include <hpp/pinocchio/device.hh>
#include <hpp/pinocchio/joint-collection.hh>
#include <hpp/pinocchio/urdf/util.hh>

Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", " << ", ";");

void test0 ()
{
  using namespace hpp::constraints;

  DevicePtr_t robot = Device::create ("robot");
  hpp::pinocchio::urdf::loadModelFromString (robot, 0, "", "anchor",
      "<robot name='foo'>"
      "  <link name='link0'>"
      "    <collision><geometry>"
      "      <sphere radius='0.1'/>"
      //"      <box size='0.1 0.1 0.1'/>"
      "    </geometry></collision>"
      "  </link>"
      "  <link name='link1'>"
      "    <collision>"
      //"      <origin rpy='0. 0.001 0.'/>"
      "      <geometry>"
      //"      <sphere radius='0.1'/>"
      //"      <sphere radius='0'/>"
      "      <box size='0.1 0.1 0.1'/>"
      "    </geometry></collision>"
      "  </link>"
      "  <joint name='joint' type='prismatic'>"
      "    <parent link='link0'/>"
      "    <child link='link1'/>"
      "    <limit effort='30' velocity='1.0' lower='-1.' upper='1.' />"
      "  </joint>"
      "</robot>", "");

  value_type r1 = 0.1, r2 = 0.05;

  ConvexCollisionAvoidancePtr_t function = ConvexCollisionAvoidance::create ("function", robot);
  value_type dub = 0.1;
  //value_type dub = 0.;
  function->addCollisionPair(0, dub, 1e-6);

  Configuration_t q (1), q_dq (1);
  value_type eps = 0.001;
  //std::cout << q << std::endl;
  //std::cout << q.size() << std::endl;

  int N = 40;
  for (int i = 0; i < N / 2; ++i)
  {
    q << ((value_type)(i+1)) / N;
    q_dq << q[0] + eps;
    value_type expectedDistance = std::min(q[0] - r1 - r2 - dub, 0.);
    value_type expectedDerivative = 2 * expectedDistance;
    expectedDistance *= expectedDistance;

    LiegroupElement f_q (function->outputSpace());
    function->value(f_q, q);
    LiegroupElement f_q_dq (function->outputSpace());
    function->value(f_q_dq, q_dq);

    matrix_t J (1, 1);
    function->jacobian(J, q);

    std::cout << q << ": " << f_q.vector() << " (" << expectedDistance << ")\n"
      << "J: " << J << " (" << (f_q_dq - f_q)/eps << ") - (" << expectedDerivative << ')' << std::endl;
    //std::cin.get();
  }
}

void test1 ()
{
  using namespace hpp::constraints;

  DevicePtr_t robot = Device::create ("robot");
  hpp::pinocchio::urdf::loadModel (robot, 0, "", "anchor",
      "/home/jmirabel/devel/hpp/src/hpp-constraints/test.urdf",
      "/home/jmirabel/devel/hpp/src/example-robot-data/robots/ur_description/srdf/ur5_gripper.srdf");

  ConvexCollisionAvoidancePtr_t function = ConvexCollisionAvoidance::create ("function", robot);
  value_type distUpperBound = 0.1,
             precision = 1e-8;

  function->addCollisionPair("shoulder_link_0", "wrist_2_link_0",
      distUpperBound, precision);
  //function->addCollisionPair("shoulder_link_0", "forearm_link_0",
      //distUpperBound, precision);

  for (int i = 0; i < 30; ++i) {
    Configuration_t q (pinocchio::randomConfiguration(robot->model()));

    LiegroupElement f (function->outputSpace());
    matrix_t J (function->outputDerivativeSize(), function->inputDerivativeSize());
    matrix_t Jfd (function->outputDerivativeSize(), function->inputDerivativeSize());
    try {
      function->value(f, q);
      function->jacobian(J, q);
      function->finiteDifferenceForward (Jfd, q, robot, 1e-4);

      if (!f.vector().isZero(1e-3)) {
        std::cout << "-------------------\n";

        std::cout << f.vector().transpose() << '\n';
        function->value(f, q);
        std::cout << f.vector().transpose() << '\n';
        //std::cout << (J - Jfd).colwise().norm().maxCoeff() << std::endl;
        std::cout << J << '\n' << Jfd << std::endl;
      }
    } catch (const std::exception& e) {
      std::cerr << "ConvexCollisionAvoidance failed: " << e.what() << '\n'
        << "q "  << q.format(CommaInitFmt) << '\n'
        << "f: " << f.vector().transpose() << '\n'
        << "J: " << J << std::endl;
    }
  }
}

int main(int argc, const char** argv)
//int main()
{
  if (argc > 1) {
    if      (strcmp(argv[1], "0") == 0) test0();
    else if (strcmp(argv[1], "1") == 0) test1();
  } else {
    test0();
    test1();
  }
  return 0;
}
