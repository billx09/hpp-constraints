//
// Copyright (c) 2016 CNRS
// Authors: Joseph Mirabel
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

#ifndef HPP_CONSTRAINTS_DIFFERENTIABLE_FUNCTION_STACK_HH
# define HPP_CONSTRAINTS_DIFFERENTIABLE_FUNCTION_STACK_HH

# include <hpp/constraints/fwd.hh>
# include <hpp/constraints/intervals.hh>
# include <hpp/constraints/differentiable-function.hh>

namespace hpp {
  namespace constraints {
    /// \addtogroup constraints
    /// \{

    /// Stack of differentiable functions
    ///
    /// This class also handles selection of cols of the output matrix.
    class HPP_CONSTRAINTS_DLLAPI DifferentiableFunctionStack :
      public DifferentiableFunction
    {
      private:
        enum CompressionType {
          VECTOR,
          MATRIX_COLUMNS,
          MATRIX
        };

      public:
        typedef std::vector<DifferentiableFunctionPtr_t> Functions_t;

        /// Return a shared pointer to a new instance
        ///
        /// \param name the name of the constraints,
        static DifferentiableFunctionStackPtr_t create (const std::string& name)
        {
          return DifferentiableFunctionStackPtr_t
            (new DifferentiableFunctionStack(name));
        }

        virtual ~DifferentiableFunctionStack () throw () {}

        /// \name Function stack management
        /// \{

        /// Get the stack of functions
        const Functions_t& functions () const
        {
          return functions_;
        }

        void add (const DifferentiableFunctionPtr_t& func)
        {
          if (functions_.empty()) {
            inputSize_           = func->inputSize();
            inputDerivativeSize_ = func->inputDerivativeSize();
            J_.resize(0, inputDerivativeSize_);
          } else {
            assert (inputSize_           == func->inputSize());
            assert (inputDerivativeSize_ == func->inputDerivativeSize());
          }
          functions_.push_back(func);
          outputSize_           += func->outputSize();
          outputDerivativeSize_ += func->outputDerivativeSize();
          J_.resize(outputDerivativeSize_, J_.cols());
        }

        /// Remove a function from the stack.
        /// \return True if the function was removed, False otherwise.
        ///
        /// \note The function stops after having removed one element. If there
        /// are duplicates, the function should be called several times.
        bool erase (const DifferentiableFunctionPtr_t& func)
        {
          for (Functions_t::iterator _func = functions_.begin();
              _func != functions_.end(); ++_func) {
            if (func == *_func) {
              functions_.erase(_func);
              outputSize_           -= func->outputSize();
              outputDerivativeSize_ -= func->outputDerivativeSize();
              J_.resize(outputDerivativeSize_, inputDerivativeSize_);
              return true;
            }
          }
          return false;
        }

        /// The output columns selection of other is not taken into account.
        void merge (const DifferentiableFunctionStackPtr_t& other)
        {
          const Functions_t& functions = other->functions();
          for (Functions_t::const_iterator _f = functions.begin();
              _f != functions.end(); ++_f)
            add (*_f);
        }

        /// \}

        /// \name Column selection management
        /// \{

        void add (const SizeInterval_t& interval)
        {
          intervals_.insert (interval);
          inputDerivativeSize_ = boost::icl::cardinality(intervals_);
        }

        SizeIntervals_t intervals () const
        {
          return intervals_;
        }

        void compressVector(vectorIn_t normal, vectorOut_t small) const
        {
          copy<true, VECTOR, vectorOut_t, vectorIn_t> (small, normal);
        }

        void uncompressVector(vectorIn_t small, vectorOut_t normal) const
        {
          copy<false, VECTOR, vectorIn_t, vectorOut_t> (small, normal);
        }

        void compressMatrix(matrixIn_t normal, matrixOut_t small, bool rows = false) const
        {
          if (rows) copy<true, MATRIX        , matrixOut_t, matrixIn_t> (small, normal);
          else      copy<true, MATRIX_COLUMNS, matrixOut_t, matrixIn_t> (small, normal);
        }

        void uncompressMatrix(matrixIn_t small, matrixOut_t normal, bool rows = false) const
        {
          if (rows) copy<false, MATRIX        , matrixIn_t, matrixOut_t> (small, normal);
          else      copy<false, MATRIX_COLUMNS, matrixIn_t, matrixOut_t> (small, normal);
        }

        /// \}

        /// Constructor
        ///
        /// \param name the name of the constraints,
        DifferentiableFunctionStack (const std::string& name)
          : DifferentiableFunction (0, 0, 0, 0, name) {}

      protected:
        void impl_compute (vectorOut_t result, ConfigurationIn_t arg) const throw ()
        {
          size_type row = 0;
          for (Functions_t::const_iterator _f = functions_.begin();
              _f != functions_.end(); ++_f) {
            const DifferentiableFunction& f = **_f;
            f.impl_compute(result.segment(row, f.outputSize()), arg);
            row += f.outputSize();
          }
        }
        void impl_jacobian (matrixOut_t jacobian, ConfigurationIn_t arg) const throw ()
        {
          bool useCache = !intervals_.empty();
          matrixOut_t J ((useCache ? matrixOut_t(J_) : jacobian));
          size_type row = 0;
          for (Functions_t::const_iterator _f = functions_.begin();
              _f != functions_.end(); ++_f) {
            const DifferentiableFunction& f = **_f;
            f.impl_jacobian(J.middleRows(row, f.outputSize()), arg);
            row += f.outputSize();
          }
          if (useCache) compressMatrix(J_, jacobian);
        }
      private:
        template <bool compress, CompressionType type, typename Small_t, typename Normal_t>
        void copy(Small_t small, Normal_t normal) const
        {
          if (intervals_.empty ()) {
            // if (compress) small = normal;
            // else          normal = small;
            assign<compress>::run(small, normal);
            return;
          }
          size_type col = 0;
          for (SizeIntervals_t::const_iterator itCol = intervals_.begin ();
              itCol != intervals_.end (); ++itCol) {
            size_type col0 = itCol->first;
            size_type nbCols = itCol->second;
            switch (type) {
              case VECTOR:
                assign<compress>::run(small .middleRows (col , nbCols), normal.middleRows (col0, nbCols));
                break;
              case MATRIX_COLUMNS:
                assign<compress>::run(small .middleCols (col , nbCols), normal.middleCols (col0, nbCols));
                break;
              case MATRIX:
                {
                  size_type row = 0;
                  for (SizeIntervals_t::const_iterator itRow = intervals_.begin ();
                      itRow != intervals_.end (); ++itRow) {
                    size_type row0 = itRow->first;
                    size_type nbRows = itRow->second;
                    assign<compress>::run(small .block (row , col , nbRows, nbCols), normal.block (row0, col0, nbRows, nbCols));
                    row += nbRows;
                  }
                  assert (row == small.rows ());
                }
                break;
            }
            col += nbCols;
          }
          assert (col == small.cols ());
        }

        template <bool rightToLeft> struct assign {
          template <typename Left_t, typename Right_t>
          static void run (Left_t l, Right_t r) { l = r; }
        };

        void resizeInput ()
        {
        }

        Functions_t functions_;
        SizeIntervals_t intervals_;

        mutable matrix_t J_;
    }; // class DifferentiableFunctionStack
    /// \}
  } // namespace constraints
} // namespace hpp

namespace hpp {
  namespace constraints {
    template <>
      struct DifferentiableFunctionStack::assign<false> {
        template <typename Left_t, typename Right_t>
          static void run (Left_t l, Right_t r) { r = l; }
      };
  } // namespace constraints
} // namespace hpp

#endif // HPP_CONSTRAINTS_DIFFERENTIABLE_FUNCTION_STACK_HH
