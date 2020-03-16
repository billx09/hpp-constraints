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

#include <hpp/constraints/convex-collision-avoidance.hh>

#include <hpp/fcl/BVH/BVH_model.h>
#include <hpp/fcl/shape/geometric_shapes.h>
#include <hpp/fcl/narrowphase/gjk.h>

#include <pinocchio/algorithm/jacobian.hpp>
#include <pinocchio/multibody/model.hpp>
#include <pinocchio/multibody/geometry.hpp>

#include <hpp/pinocchio/device.hh>
#include <hpp/pinocchio/joint-collection.hh>

namespace hpp {
namespace constraints {

ConvexCollisionAvoidance::Data::Data(value_type precision)
  : gjk (128, precision),
  epa(new EPA(128, 64, 255, precision))
{}

ConvexCollisionAvoidancePtr_t ConvexCollisionAvoidance::create (
    const std::string& name, const DevicePtr_t& robot)
{
  return ConvexCollisionAvoidancePtr_t(new ConvexCollisionAvoidance
      (name, robot));
}

ConvexCollisionAvoidance::ConvexCollisionAvoidance (const std::string& name,
          const DevicePtr_t& robot)
  : DifferentiableFunction(robot->configSize(), robot->numberDof(), 0, name),
  robot_ (robot)
{}

hpp::fcl::ShapeBase const* toShape(hpp::fcl::CollisionGeometryPtr_t geom)
{
  hpp::fcl::ShapeBase const* shape =
    dynamic_cast<hpp::fcl::ShapeBase const*> (geom.get());
  if (shape == NULL) {
    hpp::fcl::BVHModelBase const* bvh =
      dynamic_cast<hpp::fcl::BVHModelBase const*> (geom.get());
    if (bvh != NULL && bvh->convex)
      shape = dynamic_cast<hpp::fcl::ShapeBase const*> (bvh->convex.get());
  }
  if (shape == NULL) throw std::logic_error("Could not cast CollisionGeometry"
      " to ShapeBase.");

  return shape;
}

void ConvexCollisionAvoidance::addCollisionPair (::pinocchio::PairIndex index,
    value_type distanceUpperBound, value_type precision)
{
  const pinocchio::GeomModel& gmodel = robot_->geomModel();
  addCollisionPair (
      true, gmodel.collisionPairs[index].first,
      true, gmodel.collisionPairs[index].second,
      distanceUpperBound, precision);
}

void ConvexCollisionAvoidance::addCollisionPair (const std::string& geoma,
    const std::string& geomb,
    value_type distanceUpperBound, value_type precision)
{
  const pinocchio::GeomModel& gmodel = robot_->geomModel();
  bool robota = true;
  bool robotb = true;
  if (!gmodel.existGeometryName(geoma)) {
    robota = false;
    if (!obsModel_->existGeometryName(geoma))
      throw std::invalid_argument("Could not find geometry " + geoma);
  }
  if (!gmodel.existGeometryName(geomb)) {
    robotb = false;
    if (!obsModel_->existGeometryName(geomb))
    throw std::invalid_argument("Could not find geometry " + geomb);
  }
  if (!robota && !robotb)
    throw std::invalid_argument("Both geometry are obstacles");

  addCollisionPair (
      robota, robota ? gmodel.getGeometryId(geoma) : obsModel_->getGeometryId(geoma),
      robotb, robotb ? gmodel.getGeometryId(geomb) : obsModel_->getGeometryId(geomb),
      distanceUpperBound, precision);
}

void ConvexCollisionAvoidance::addCollisionPair (::pinocchio::GeomIndex ida,
    ::pinocchio::GeomIndex idb,
    value_type distanceUpperBound, value_type precision)
{
  addCollisionPair (true, ida, true, idb, distanceUpperBound, precision);
}

void ConvexCollisionAvoidance::addCollisionPair (
    bool robota, ::pinocchio::GeomIndex ida,
    bool robotb, ::pinocchio::GeomIndex idb,
    value_type distanceUpperBound, value_type precision)
{
  const pinocchio::GeomModel& gmodel = robot_->geomModel();

  fcl::Transform3f tfa, tfb;
  tfa.setIdentity();
  tfb.setIdentity();
  tfb.setTranslation(vector3_t::Ones());

  Data data (precision);
  data.ida = ida;
  data.aRob = robota;
  data.idb = idb;
  data.bRob = robotb;

  const ::pinocchio::GeometryObject
    ga = (robota ? gmodel : *obsModel_).geometryObjects[data.ida],
    gb = (robotb ? gmodel : *obsModel_).geometryObjects[data.idb];

  if (robota)
    std::cout << ga.name << ": " << robot_->model().names[ga.parentJoint] << std::endl;

  data.shape.set (toShape(ga.geometry), toShape(gb.geometry), tfa, tfb);

  // Stop if you can prove they are more than 0.1m away from each other.
  data.gjk.setDistanceEarlyBreak(distanceUpperBound);
  data.gjk.ray << 1, 0, 0;
  data.distanceUpperBound = distanceUpperBound;

  datas_.push_back(data);

  // Update the function output space
  outputSpace_ = LiegroupSpace::Rn(datas_.size());
}

void ConvexCollisionAvoidance::impl_compute (LiegroupElementRef result,
    ConfigurationIn_t argument) const
{
  computeCollision(argument);

  value_type dist;
  for (std::size_t i = 0; i < datas_.size(); ++i) {
    dist = datas_[i].distance - datas_[i].distanceUpperBound;
    //result.vector()[i] = (dist < 0 ? dist * dist : 0);
    result.vector()[i] = std::min(dist, 0.);
  }
}

void ConvexCollisionAvoidance::impl_jacobian (matrixOut_t jacobian,
    ConfigurationIn_t arg) const
{
  computeCollision(arg);

  hpp::pinocchio::DeviceSync device (robot_);
  device.currentConfiguration (arg);
  device.computeForwardKinematics ();
  device.updateGeometryPlacements ();

  value_type dist;
  const pinocchio::GeomModel& gmodel = device.geomModel();
  for (std::size_t i = 0; i < datas_.size(); ++i) {
    Data& data = datas_[i];

    dist = data.distance - data.distanceUpperBound;
    if (dist >= 0) {
      jacobian.row(i).setZero();
      continue;
    }

    const ::pinocchio::GeometryObject
      ga = (data.aRob ? gmodel : *obsModel_).geometryObjects[data.ida],
      gb = (data.bRob ? gmodel : *obsModel_).geometryObjects[data.idb];

    // dist = || w0 - w1 ||^2
    // w = w0 - w1
    // J = 2 * w^T * 0J01
    // a = j0, b = j1

    //Transform3f oM1 (data.shape.oR1, data.shape.ot1);
    static const Transform3f Id (Transform3f::Identity());
    Transform3f aMb (
        (data.aRob ? device.data().oMi[ga.parentJoint] : Id).actInv(
        (data.bRob ? device.data().oMi[gb.parentJoint] : Id)));

    vector3_t normal = ga.placement.rotation() * data.normal;
    //std::cout << data.normal.transpose() << std::endl;
    //std::cout <<      normal.transpose() << std::endl;

    //JointJacobian_t aJa(6, robot_->numberDof()), bJb(6, robot_->numberDof()),
                    //J(6, robot_->numberDof());
    JointJacobian_t J(6, robot_->numberDof()), bJb(6, robot_->numberDof());
    J.setZero();

    if (data.aRob && ga.parentJoint > 0)
      ::pinocchio::getJointJacobian (device.model(), device.data(), ga.parentJoint,
          ::pinocchio::LOCAL, J);
    if (data.bRob && gb.parentJoint > 0) {
      bJb.setZero();
      ::pinocchio::getJointJacobian (device.model(), device.data(), gb.parentJoint,
          ::pinocchio::LOCAL, bJb);
      J -= aMb.toActionMatrix() * bJb;
    }

    // J = aJa - aXb bJb
    //JointJacobian_t J (aJa - aMb.toActionMatrix() * bJb);
    // JointJacobian_t J (6, aJa.cols());
    // ::pinocchio::motionSet(aMb, bJb, J);
    // J *= -1;
    // J += aJa;

    //jacobian.row(i) = - (2 * (data.distance - data.distanceUpperBound))
      //* data.normal.transpose() * (ga.placement.toActionMatrixInverse() * J).topRows<3>();
    jacobian.row(i) = - data.normal.transpose() * (ga.placement.toActionMatrixInverse() * J).topRows<3>();
  }
}

void ConvexCollisionAvoidance::computeCollision(ConfigurationIn_t arg) const
{
  hpp::pinocchio::DeviceSync device (robot_);
  device.currentConfiguration (arg);
  device.computeForwardKinematics ();
  device.updateGeometryPlacements ();

  //const pinocchio::GeomModel& gmodel = device.geomModel();
  const pinocchio::GeomData & gdata  = device.geomData ();
  for (std::size_t i = 0; i < datas_.size(); ++i) {
    Data& data = datas_[i];

    // const ::pinocchio::GeometryObject
    //   ga = gmodel.geometryObjects[data.ida],
    //   gb = gmodel.geometryObjects[data.idb];

    // Update
    Transform3f oM1 (
        (data.aRob ? gdata : *obsData_).oMg[data.ida].actInv(
        (data.bRob ? gdata : *obsData_).oMg[data.idb]));
    data.shape.oR1 = oM1.rotation();
    data.shape.ot1 = oM1.translation();

    //data.shape.set (toShape(ga.geometry), toShape(gb.geometry),
    //    ::pinocchio::toFclTransform3f(gdata.oMg[data.ida]),
    //    ::pinocchio::toFclTransform3f(gdata.oMg[data.idb]));

    // Call GJK
    data.gjk_status = data.gjk.evaluate(
        data.shape, data.gjk.getGuessFromSimplex());

    switch(data.gjk_status) {
      case GJK::Inside:
        if (data.gjk.distance == 0 && !data.gjk.ray.isZero()) {
          // Distance is very low but object are not colliding.
          // Don't run EPA
          data.distance = data.gjk.distance;
          data.gjk.getClosestPoints (data.shape, data.w0, data.w1);
          // TODO The GJK ray should provide a more robust normal.
          data.normal = (data.w1 - data.w0).normalized();
          break;
        }
        //TODO ray == 0 at this stage.
        data.epa_status = data.epa->evaluate(data.gjk, data.gjk.ray);
        switch (data.epa_status) {
          case EPA::OutOfFaces:
          case EPA::OutOfVertices:
            throw std::runtime_error("EPA ran out of memory.");
          case EPA::Failed:
            throw std::runtime_error("EPA failed: unknown error.");
          case EPA::InvalidHull:
            throw std::runtime_error("EPA found an invalid hull.");
          case EPA::NonConvex:
            throw std::runtime_error("EPA ran on a non-convex model.");
          case EPA::AccuracyReached:
          case EPA::Valid:
            data.distance = - data.epa->depth;
            data.normal = data.epa->normal;
            data.epa->getClosestPoints (data.shape, data.w0, data.w1);
            break;
          //enum Status {Touching, Degenerated, FallBack};
          default:
            throw std::runtime_error("EPA failed: unhandled error code.");
        }
        break;
      case GJK::Valid:
        data.distance = data.gjk.distance;
        if (data.distance < data.distanceUpperBound) {
          data.gjk.getClosestPoints (data.shape, data.w0, data.w1);
          // TODO The GJK ray should provide a more robust normal.
          data.normal = (data.w1 - data.w0) / data.gjk.distance;
        }
        break;
      case GJK::Failed:
        throw std::runtime_error("GJK failed.");
    }
  }
}

} // namespace constraints
} // namespace hpp
