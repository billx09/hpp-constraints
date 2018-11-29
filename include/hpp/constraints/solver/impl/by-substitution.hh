// Copyright (c) 2017, Joseph Mirabel
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

#ifndef HPP_CONSTRAINTS_SOLVER_IMPL_BY_SUBSTITUTION_HH
#define HPP_CONSTRAINTS_SOLVER_IMPL_BY_SUBSTITUTION_HH

namespace hpp {
  namespace constraints {
    namespace solver {
    template <typename LineSearchType>
    inline HierarchicalIterative::Status BySubstitution::impl_solve (
        vectorOut_t arg,
        LineSearchType lineSearch) const
    {
      assert (!arg.hasNaN());

      explicit_.solve(arg);

      size_type errorDecreased = 3, iter = 0;
      value_type previousSquaredNorm =
	std::numeric_limits<value_type>::infinity();
      static const value_type dqMinSquaredNorm = Eigen::NumTraits<value_type>::dummy_precision();
      value_type initSquaredNorm = 0;

      // Fill value and Jacobian
      computeValue<true> (arg);
      computeError();

      bool errorWasBelowThr = (squaredNorm_ < squaredErrorThreshold_);
      vector_t initArg;
      if (errorWasBelowThr) {
        initArg = arg;
        iter = std::max (maxIterations_,size_type(2)) - 2;
        initSquaredNorm = squaredNorm_;
      }

      if (squaredNorm_ > .25 * squaredErrorThreshold_
          && reducedDimension_ == 0) return INFEASIBLE;

      Status status;
      while (squaredNorm_ > .25 * squaredErrorThreshold_ && errorDecreased &&
	     iter < maxIterations_) {

        // Update the jacobian using the jacobian of the explicit system.
        updateJacobian(arg);
        computeSaturation(arg);
        computeDescentDirection ();
        if (dq_.squaredNorm () < dqMinSquaredNorm) {
          // TODO INFEASIBLE means that we have reached a local minima.
          // The problem may still be feasible from a different starting point.
          status = INFEASIBLE;
          break;
        }
        lineSearch (*this, arg, dq_);
        explicit_.solve(arg);

	computeValue<true> (arg);
        computeError ();

	--errorDecreased;
	if (squaredNorm_ < previousSquaredNorm)
          errorDecreased = 3;
        else
          status = ERROR_INCREASED;
	previousSquaredNorm = squaredNorm_;
	++iter;

      }

      if (errorWasBelowThr) {
        if (squaredNorm_ > initSquaredNorm) {
          arg = initArg;
        }
        return SUCCESS;
      }

      if (squaredNorm_ > squaredErrorThreshold_) {
        return (iter >= maxIterations_) ? MAX_ITERATION_REACHED : status;
      }
      assert (!arg.hasNaN());
      return SUCCESS;
    }

    template <typename LineSearchType>
    inline HierarchicalIterative::Status BySubstitution::impl_solve (
        vectorOut_t arg,
        vectorOut_t v,
        LineSearchType lineSearch) const
    {
      assert (!arg.hasNaN());
      if (datas_.size() != 1) {
        throw std::runtime_error ("Projection of a velocity using a BySubstitution"
            " only works with 1 level");
      }

      explicit_.solve(arg);

      size_type errorDecreased = 3, iter = 0;
      value_type previousSquaredNorm =
	std::numeric_limits<value_type>::infinity();
      static const value_type dqMinSquaredNorm = Eigen::NumTraits<value_type>::dummy_precision();
      value_type initSquaredNorm = 0;

      // Fill value and Jacobian
      computeValue<true> (arg);
      computeError();

      bool errorWasBelowThr = (squaredNorm_ < squaredErrorThreshold_);
      vector_t initArg;
      if (errorWasBelowThr) {
        initArg = arg;
        iter = std::max (maxIterations_,size_type(2)) - 2;
        initSquaredNorm = squaredNorm_;
      }

      if (squaredNorm_ > .25 * squaredErrorThreshold_
          && reducedDimension_ == 0) return INFEASIBLE;

      const Data& d = datas_[0];
      vector_t vn (freeVariables_.rview(v)),
               J_times_vn (d.reducedJ.rows());

      Status status;
      while (squaredNorm_ > .25 * squaredErrorThreshold_ && errorDecreased &&
	     iter < maxIterations_) {

        // Update the jacobian using the jacobian of the explicit system.
        updateJacobian(arg);
        computeSaturation(arg);
        computeDescentDirection ();
        if (dq_.squaredNorm () < dqMinSquaredNorm) {
          // TODO INFEASIBLE means that we have reached a local minima.
          // The problem may still be feasible from a different starting point.
          status = INFEASIBLE;
          break;
        }
        value_type alpha = lineSearch (*this, arg, dq_);
        explicit_.solve(arg);
        // Compute v_n
        J_times_vn = d.reducedJ * vn;
        vn -= alpha * datas_[0].svd.solve(J_times_vn);

	computeValue<true> (arg);
        computeError ();

	--errorDecreased;
	if (squaredNorm_ < previousSquaredNorm)
          errorDecreased = 3;
        else
          status = ERROR_INCREASED;
	previousSquaredNorm = squaredNorm_;
	++iter;

      }

      freeVariables_.lview(v) = vn;
      explicit_.inDers().transpose().lview(v) =
        ExplicitConstraintSet::MatrixBlockView(JeExpanded_,
            explicit_.inDers().nbIndices(), explicit_.inDers().indices(),
            explicit_.inDers().nbIndices(), explicit_.inDers().indices()
            ).eval() * explicit_.inDers().transpose().rview(v).eval();

      if (errorWasBelowThr) {
        if (squaredNorm_ > initSquaredNorm) {
          arg = initArg;
        }
        return SUCCESS;
      }

      if (squaredNorm_ > squaredErrorThreshold_) {
        return (iter >= maxIterations_) ? MAX_ITERATION_REACHED : status;
      }
      assert (!arg.hasNaN());
      return SUCCESS;
    }
    } // namespace solver
  } // namespace constraints
} // namespace hpp

#endif // HPP_CONSTRAINTS_SOLVER_IMPL_BY_SUBSTITUTION_HH
