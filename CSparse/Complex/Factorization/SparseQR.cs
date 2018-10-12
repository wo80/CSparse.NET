// -----------------------------------------------------------------------
// <copyright file="SparseQR.cs">
// Copyright (c) 2006-2016, Timothy A. Davis
// Copyright (c) 2012-2016, Christian Woltering
// </copyright>
// -----------------------------------------------------------------------

namespace CSparse.Complex.Factorization
{
    using CSparse.Factorization;
    using CSparse.Ordering;
    using CSparse.Storage;
    using System;
    using System.Numerics;

    /// <summary>
    /// Sparse QR decomposition.
    /// </summary>
    /// <remarks>
    /// See Chapter 5 (Orthogonal methods) in "Direct Methods for Sparse Linear Systems"
    /// by Tim Davis.
    /// </remarks>
    public class SparseQR : SparseQR<Complex>
    {
        #region Static methods

        /// <summary>
        /// Creates a sparse QR factorization.
        /// </summary>
        /// <param name="A">Column-compressed matrix, symmetric positive definite.</param>
        /// <param name="order">Ordering method to use (natural or A+A').</param>
        public static SparseQR Create(CompressedColumnStorage<Complex> A, ColumnOrdering order)
        {
            return Create(A, order, null);
        }

        /// <summary>
        /// Creates a sparse QR factorization.
        /// </summary>
        /// <param name="A">Column-compressed matrix, symmetric positive definite.</param>
        /// <param name="order">Ordering method to use (natural or A+A').</param>
        /// <param name="progress">Report progress (range from 0.0 to 1.0).</param>
        public static SparseQR Create(CompressedColumnStorage<Complex> A, ColumnOrdering order,
            IProgress<double> progress)
        {
            Check.NotNull(A, "A");

            int m = A.RowCount;
            int n = A.ColumnCount;

            var C = new SparseQR(m, n);

            if (m >= n)
            {
                var p = AMD.Generate(A, order);

                // Ordering and symbolic analysis
                C.SymbolicAnalysis(A, p, order == ColumnOrdering.Natural);

                // Numeric QR factorization
                C.Factorize(A, progress);
            }
            else
            {
                // Ax=b is underdetermined
                var AT = A.Transpose();

                var p = AMD.Generate(AT, order);

                // Ordering and symbolic analysis
                C.SymbolicAnalysis(AT, p, order == ColumnOrdering.Natural);

                // Numeric QR factorization of A'
                C.Factorize(AT, progress);
            }

            return C;
        }

        #endregion
        
        private SparseQR(int rows, int columns)
            : base(rows, columns)
        {
        }

        /// <summary>
        /// Solves a system of linear equations, <c>Ax = b</c>.
        /// </summary>
        /// <param name="input">The right hand side vector, <c>b</c>.</param>
        /// <param name="result">The left hand side vector, <c>x</c>.</param>
        /// <remarks>
        /// Let A be a m-by-n matrix. If m >= n a least-squares problem (min |Ax-b|)
        /// is solved. If m &lt; n the underdetermined system is solved.
        /// </remarks>
        public override void Solve(Complex[] input, Complex[] result)
        {
            if (input == null) throw new ArgumentNullException("input");

            if (result == null) throw new ArgumentNullException("result");

            var x = new Complex[S.m2];

            if (m >= n)
            {
                // x(0:m-1) = b(p(0:m-1)
                Permutation.ApplyInverse(S.pinv, input, x, m);

                // Apply Householder reflection to x.
                for (int k = 0; k < n; k++)
                {
                    ApplyHouseholder(Q, k, beta[k], x);
                }

                SolverHelper.SolveUpper(R, x); // x = R\x

                // b(q(0:n-1)) = x(0:n-1)
                Permutation.ApplyInverse(S.q, x, result, n);
            }
            else
            {
                // x(q(0:m-1)) = b(0:m-1)
                Permutation.Apply(S.q, input, x, m);

                SolverHelper.SolveUpperTranspose(R, x); // x = R'\x

                // Apply Householder reflection to x.
                for (int k = m - 1; k >= 0; k--)
                {
                    ApplyHouseholder(Q, k, beta[k], x);
                }

                // b(0:n-1) = x(p(0:n-1))
                Permutation.Apply(S.pinv, x, result, n);
            }
        }

        /// <summary>
        /// Create a Householder reflection [v,beta,s]=house(x), overwrite x with v,
        /// where (I-beta*v*v')*x = s*e1 and e1 = [1 0 ... 0]'.
        /// </summary>
        /// <remarks>
        /// Note that this CXSparse version is different than CSparse.  See Higham,
        /// Accuracy and Stability of Num Algorithms, 2nd ed, 2002, page 357.
        /// </remarks>
        protected override Complex CreateHouseholder(Complex[] x, int offset, ref double beta, int n)
        {
            Complex s = Complex.Zero;
            int i;
            if (x == null) return -1; // check inputs

            // s = norm(x)
            for (i = 0; i < n; i++)
            {
                s += x[offset + i] * Complex.Conjugate(x[offset + i]);
            }

            s = Complex.Sqrt(s);
            if (s == 0)
            {
                beta = 0;
                x[offset] = 1;
            }
            else
            {
                // s = sign(x[0]) * norm (x) ;
                if (x[offset] != 0)
                {
                    s *= x[offset] / Complex.Abs(x[offset]);
                }
                x[offset] += s;
                beta = 1 / (Complex.Conjugate(s) * x[offset]).Real;
            }
            return (-s);
        }

        /// <summary>
        /// Apply the ith Householder vector to x.
        /// </summary>
        protected override bool ApplyHouseholder(CompressedColumnStorage<Complex> V, int i, double beta, Complex[] x)
        {
            int p = 0;
            Complex tau = Complex.Zero;

            if (x == null) return false; // check inputs

            var vp = V.ColumnPointers;
            var vi = V.RowIndices;
            var vx = V.Values;

            var vpi1 = vp[i + 1];

            for (p = vp[i]; p < vpi1; p++) // tau = v'*x
            {
                tau += Complex.Conjugate(vx[p]) * x[vi[p]];
            }
            tau *= beta; // tau = beta*(v'*x)
            for (p = vp[i]; p < vpi1; p++) // x = x - v*tau
            {
                x[vi[p]] -= vx[p] * tau;
            }
            return true;
        }
    }
}
