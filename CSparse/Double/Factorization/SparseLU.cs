namespace CSparse.Double.Factorization
{
    using CSparse.Factorization;
    using CSparse.Ordering;
    using CSparse.Properties;
    using CSparse.Storage;
    using System;

    /// <summary>
    /// Sparse LU decomposition.
    /// </summary>
    /// <remarks>
    /// See Chapter 6 (LU factorization) in "Direct Methods for Sparse Linear Systems"
    /// by Tim Davis.
    /// </remarks>
    public class SparseLU : ISparseFactorization<double>
    {
        readonly int n;

        SymbolicFactorization S;
        CompressedColumnStorage<double> L, U;
        int[] pinv; // partial pivoting

        double[] temp; // workspace
        
        #region Static methods

        /// <summary>
        /// Creates a LU factorization.
        /// </summary>
        /// <param name="A">Column-compressed matrix, symmetric positive definite.</param>
        /// <param name="order">Ordering method to use (natural or A+A').</param>
        /// <param name="tol">Partial pivoting tolerance (form 0.0 to 1.0).</param>
        public static SparseLU Create(CompressedColumnStorage<double> A, ColumnOrdering order,
            double tol)
        {
            return Create(A, order, tol, null);
        }

        /// <summary>
        /// Creates a LU factorization.
        /// </summary>
        /// <param name="A">Column-compressed matrix, symmetric positive definite.</param>
        /// <param name="order">Ordering method to use (natural or A+A').</param>
        /// <param name="tol">Partial pivoting tolerance (form 0.0 to 1.0).</param>
        /// <param name="progress">Report progress (range from 0.0 to 1.0).</param>
        public static SparseLU Create(CompressedColumnStorage<double> A, ColumnOrdering order,
            double tol, IProgress<double> progress)
        {
            return Create(A, AMD.Generate(A, order), tol, progress);
        }

        /// <summary>
        /// Creates a LU factorization.
        /// </summary>
        /// <param name="A">Column-compressed matrix, symmetric positive definite.</param>
        /// <param name="p">Permutation.</param>
        /// <param name="tol">Partial pivoting tolerance (form 0.0 to 1.0).</param>
        public static SparseLU Create(CompressedColumnStorage<double> A, int[] p, double tol)
        {
            return Create(A, p, tol, null);
        }

        /// <summary>
        /// Creates a LU factorization.
        /// </summary>
        /// <param name="A">Column-compressed matrix, symmetric positive definite.</param>
        /// <param name="p">Permutation.</param>
        /// <param name="tol">Partial pivoting tolerance (form 0.0 to 1.0).</param>
        /// <param name="progress">Report progress (range from 0.0 to 1.0).</param>
        public static SparseLU Create(CompressedColumnStorage<double> A, int[] p, double tol,
            IProgress<double> progress)
        {
            Check.NotNull(A, "A");
            Check.NotNull(p, "p");

            int n = A.ColumnCount;

            Check.SquareMatrix(A, "A");
            Check.Permutation(p, n, "p");
            Check.NotNaN(tol, "tol");

            // Ensure tol is in range.
            tol = Math.Min(Math.Max(tol, 0.0), 1.0);

            var C = new SparseLU(n);

            // Ordering and symbolic analysis
            C.SymbolicAnalysis(A, p);

            // Numeric Cholesky factorization
            C.Factorize(A, tol, progress);

            return C;
        }

        #endregion

        private SparseLU(int n)
        {
            this.n = n;
            this.temp = new double[n];
        }

        /// <summary>
        /// Gets the number of nonzeros in both L and U factors together.
        /// </summary>
        public int NonZerosCount
        {
            get { return (L.NonZerosCount + U.NonZerosCount - n); }
        }

        /// <summary>
        /// Solves a system of linear equations, <c>Ax = b</c>.
        /// </summary>
        /// <param name="input">The right hand side vector, <c>b</c>.</param>
        /// <param name="result">The left hand side vector, <c>x</c>.</param>
        public void Solve(double[] input, double[] result)
        {
            if (input == null) throw new ArgumentNullException(nameof(input));

            if (result == null) throw new ArgumentNullException(nameof(result));

            var x = this.temp;

            Permutation.ApplyInverse(pinv, input, x, n); // x = b(p)

            SolverHelper.SolveLower(L, x); // x = L\x.

            SolverHelper.SolveUpper(U, x); // x = U\x.

            Permutation.ApplyInverse(S.q, x, result, n); // b(q) = x
        }

        /// <summary>
        /// Solves a system of linear equations, <c>A'x = b</c>.
        /// </summary>
        /// <param name="input">The right hand side vector, <c>b</c>.</param>
        /// <param name="result">The left hand side vector, <c>x</c>.</param>
        public void SolveTranspose(double[] input, double[] result)
        {
            if (input == null) throw new ArgumentNullException(nameof(input));

            if (result == null) throw new ArgumentNullException(nameof(result));

            var x = this.temp;

            Permutation.Apply(S.q, input, x, n); // x = Q'*b

            SolverHelper.SolveUpperTranspose(U, x); // x = U'\x.

            SolverHelper.SolveLowerTranspose(L, x); // x = L'\x.

            Permutation.Apply(pinv, x, result, n); // b = P'*x
        }

        /// <summary>
        /// [L,U,pinv] = lu(A, [q lnz unz]). lnz and unz can be guess.
        /// </summary>
        private void Factorize(CompressedColumnStorage<double> A, double tol, IProgress<double> progress)
        {
            int[] q = S.q;

            int i;
            int lnz = S.lnz;
            int unz = S.unz;

            this.L = CompressedColumnStorage<double>.Create(n, n, lnz);
            this.U = CompressedColumnStorage<double>.Create(n, n, unz);
            this.pinv = new int[n];

            // Workspace
            var x = this.temp;
            var xi = new int[2 * n];

            for (i = 0; i < n; i++)
            {
                // No rows pivotal yet.
                pinv[i] = -1;
            }

            lnz = unz = 0;

            int ipiv, top, p, col;
            double pivot;
            double a, t;

            int[] li, ui;
            int[] lp = L.ColumnPointers;
            int[] up = U.ColumnPointers;
            double[] lx, ux;

            double current = 0.0;
            double step = n / 100.0;

            // Now compute L(:,k) and U(:,k)
            for (int k = 0; k < n; k++)
            {
                // Progress reporting.
                if (k >= current)
                {
                    current += step;

                    if (progress != null)
                    {
                        progress.Report(k / (double)n);
                    }
                }

                // Triangular solve
                lp[k] = lnz; // L(:,k) starts here
                up[k] = unz; // U(:,k) starts here

                if (lnz + n > L.Values.Length) L.Resize(2 * L.Values.Length + n);
                if (unz + n > U.Values.Length) U.Resize(2 * U.Values.Length + n);

                li = L.RowIndices;
                ui = U.RowIndices;
                lx = L.Values;
                ux = U.Values;
                col = q != null ? (q[k]) : k;
                top = SolveSp(L, A, col, xi, x, pinv, true);  // x = L\A(:,col)

                // Find pivot
                ipiv = -1;
                a = -1;
                for (p = top; p < n; p++)
                {
                    i = xi[p]; // x(i) is nonzero
                    if (pinv[i] < 0) // Row i is not yet pivotal
                    {
                        if ((t = Math.Abs(x[i])) > a)
                        {
                            a = t; // Largest pivot candidate so far
                            ipiv = i;
                        }
                    }
                    else // x(i) is the entry U(pinv[i],k)
                    {
                        ui[unz] = pinv[i];
                        ux[unz++] = x[i];
                    }
                }

                if (ipiv == -1 || a <= 0.0)
                {
                    throw new Exception("No pivot element found.");
                }

                if (pinv[col] < 0 && Math.Abs(x[col]) >= a * tol)
                {
                    ipiv = col;
                }

                // Divide by pivot
                pivot = x[ipiv]; // the chosen pivot
                ui[unz] = k; // last entry in U(:,k) is U(k,k)
                ux[unz++] = pivot;
                pinv[ipiv] = k; // ipiv is the kth pivot row
                li[lnz] = ipiv; // first entry in L(:,k) is L(k,k) = 1
                lx[lnz++] = 1.0;
                for (p = top; p < n; p++) // L(k+1:n,k) = x / pivot
                {
                    i = xi[p];
                    if (pinv[i] < 0) // x(i) is an entry in L(:,k)
                    {
                        li[lnz] = i; // save unpermuted row in L
                        lx[lnz++] = x[i] / pivot; // scale pivot column
                    }
                    x[i] = 0.0; // x [0..n-1] = 0 for next k
                }
            }

            // Finalize L and U
            lp[n] = lnz;
            up[n] = unz;
            li = L.RowIndices; // fix row indices of L for final pinv
            for (p = 0; p < lnz; p++)
            {
                li[p] = pinv[li[p]];
            }

            // Remove extra space from L and U
            L.Resize(0);
            U.Resize(0);
        }

        /// <summary>
        /// Symbolic ordering and analysis for LU.
        /// </summary>
        /// <param name="A"></param>
        /// <param name="p">Permutation.</param>
        private void SymbolicAnalysis(CompressedColumnStorage<double> A, int[] p)
        {
            var sym = this.S = new SymbolicFactorization();

            // Fill-reducing ordering
            sym.q = p;
            
            // Guess nnz(L) and nnz(U)
            sym.unz = sym.lnz = 4 * (A.ColumnPointers[n]) + n;
        }

        /// <summary>
        /// Solve Gx=b(:,k), where G is either upper (lo=false) or lower (lo=true)
        /// triangular.
        /// </summary>
        /// <param name="G">lower or upper triangular matrix in column-compressed form</param>
        /// <param name="B">right hand side, b=B(:,k)</param>
        /// <param name="k">use kth column of B as right hand side</param>
        /// <param name="xi">size 2*n, nonzero pattern of x in xi[top..n-1]</param>
        /// <param name="x">size n, x in x[xi[top..n-1]]</param>
        /// <param name="pinv">mapping of rows to columns of G, ignored if null</param>
        /// <param name="lo">true if lower triangular, false if upper</param>
        /// <returns>top, -1 in error</returns>
        private int SolveSp(CompressedColumnStorage<double> G, CompressedColumnStorage<double> B,
            int k, int[] xi, double[] x, int[] pinv, bool lo)
        {
            if (xi == null || x == null) return -1;

            var gp = G.ColumnPointers;
            var gi = G.RowIndices;
            var gx = G.Values;

            var bp = B.ColumnPointers;
            var bi = B.RowIndices;
            var bx = B.Values;

            int n = G.ColumnCount;

            // xi[top..n-1]=Reach(B(:,k))
            int top = GraphHelper.Reach(gp, gi, bp, bi, n, k, xi, pinv);

            int j, J, p, q, px;

            for (p = top; p < n; p++)
            {
                x[xi[p]] = 0; // clear x
            }

            for (p = bp[k]; p < bp[k + 1]; p++)
            {
                x[bi[p]] = bx[p]; // scatter B
            }

            for (px = top; px < n; px++)
            {
                j = xi[px]; // x(j) is nonzero
                J = pinv != null ? (pinv[j]) : j; // j maps to col J of G
                if (J < 0) continue; // column J is empty
                x[j] /= gx[lo ? (gp[J]) : (gp[J + 1] - 1)]; // x(j) /= G(j,j)
                p = lo ? (gp[J] + 1) : (gp[J]); // lo: L(j,j) 1st entry
                q = lo ? (gp[J + 1]) : (gp[J + 1] - 1); // up: U(j,j) last entry
                for (; p < q; p++)
                {
                    x[gi[p]] -= gx[p] * x[j]; // x(i) -= G(i,j) * x(j)
                }
            }

            // Return top of stack.
            return top;
        }
    }
}
