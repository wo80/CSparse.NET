// -----------------------------------------------------------------------
// <copyright file="SparseLDL.cs">
// Copyright (c) 2005-2013, Timothy A. Davis
// Copyright (c) 2012-2016, Christian Woltering
// </copyright>
// -----------------------------------------------------------------------

namespace CSparse.Double.Factorization
{
    using CSparse.Factorization;
    using CSparse.Ordering;
    using CSparse.Storage;
    using System;

    /// <summary>
    /// Sparse LDL' factorization.
    /// </summary>
    /// <remarks>
    /// If A is positive definite then the factorization will be accurate. A can be
    /// indefinite (with negative values on the diagonal D), but in this case no
    /// guarantee of accuracy is provided, since no numeric pivoting is performed.
    ///
    /// Only the diagonal and upper triangular part of A (or PAP' if a permutation
    /// P is provided) is accessed.  The lower triangular parts of the matrix A or
    /// PAP' can be present, but they are ignored.
    /// </remarks>
    public class SparseLDL : ISparseFactorization<double>
    {
        SymbolicFactorization S;
        CompressedColumnStorage<double> L;
        double[] D;

        int n;

        public SparseLDL(CompressedColumnStorage<double> A, ColumnOrdering order)
        {
            if ((int)order > 1) // AtA ordering not allowed
            {
                throw new ArgumentException("order");
            }

            this.n = A.ColumnCount;

            // Ordering and symbolic analysis
            SymbolicAnalysis(order, A);
            
            // Numeric Cholesky factorization
            Factorize(A);
        }

        /// <summary>
        /// Gets the number of nonzeros of the L.
        /// </summary>
        public int NonZerosCount
        {
            get { return L.NonZerosCount; }
        }

        /// <summary>
        /// Solves a linear system Ax=b, where A is symmetric positive definite.
        /// </summary>
        /// <param name="input">Right hand side b.</param>
        /// <param name="result">Solution vector x.</param>
        public void Solve(double[] input, double[] result)
        {
            if (input == null) throw new ArgumentNullException("input");

            if (result == null) throw new ArgumentNullException("result");

            double[] x = new double[n];

            Permutation.ApplyInverse(S.pinv, input, x, n); // x = P*b

            var lx = L.Values;
            var lp = L.ColumnPointers;
            var li = L.RowIndices;

            var d = this.D;

            int end;

            // Solve lower triangular system by forward elimination, x = L\x.
            for (int i = 0; i < n; i++)
            {
                end = lp[i + 1];
                for (int p = lp[i]; p < end; p++)
                {
                    x[li[p]] -= lx[p] * x[i];
                }
            }

            // Solve diagonal system, x = D\x.
            for (int i = 0; i < n; i++)
            {
                x[i] /= d[i];
            }

            // Solve upper triangular system by backward elimination, x = L'\x.
            for (int i = n - 1; i >= 0; i--)
            {
                end = lp[i + 1];
                for (int p = lp[i]; p < end; p++)
                {
                    x[i] -= lx[p] * x[li[p]];
                }
            }

            Permutation.Apply(S.pinv, x, result, n); // b = P'*x
        }

        /// <summary>
        /// Ordering and symbolic analysis for a LDL' factorization.
        /// </summary>
        /// <param name="order">Column ordering.</param>
        /// <param name="A">Matrix to factorize.</param>
        private void SymbolicAnalysis(ColumnOrdering order, CompressedColumnStorage<double> A)
        {
            int n = A.ColumnCount;

            var sym = this.S = new SymbolicFactorization();

            var ap = A.ColumnPointers;
            var ai = A.RowIndices;

            // P = amd(A+A') or natural
            var P = AMD.Generate(A, order);
            var Pinv = Permutation.Invert(P);

            // Output: column pointers and elimination tree.
            var lp = new int[n + 1];
            var parent = new int[n];

            // Workspace
            var lnz = new int[n];
            var flag = new int[n];

            int i, k, p, kk, p2;

            for (k = 0; k < n; k++)
            {
                // L(k,:) pattern: all nodes reachable in etree from nz in A(0:k-1,k) 
                parent[k] = -1; // parent of k is not yet known 
                flag[k] = k; // mark node k as visited 
                lnz[k] = 0; // count of nonzeros in column k of L 
                kk = (P != null) ? (P[k]) : (k); // kth original, or permuted, column 
                p2 = ap[kk + 1];
                for (p = ap[kk]; p < p2; p++)
                {
                    // A(i,k) is nonzero (original or permuted A) 
                    i = (Pinv != null) ? (Pinv[ai[p]]) : (ai[p]);
                    if (i < k)
                    {
                        // follow path from i to root of etree, stop at flagged node 
                        for (; flag[i] != k; i = parent[i])
                        {
                            // find parent of i if not yet determined 
                            if (parent[i] == -1) parent[i] = k;
                            lnz[i]++; // L(k,i) is nonzero 
                            flag[i] = k; // mark i as visited 
                        }
                    }
                }
            }

            // construct Lp index array from Lnz column counts 
            lp[0] = 0;
            for (k = 0; k < n; k++)
            {
                lp[k + 1] = lp[k] + lnz[k];
            }

            sym.parent = parent;
            sym.cp = lp;
            sym.q = P;
            sym.pinv = Pinv;
        }

        /// <summary>
        /// Compute the numeric LDL' factorization of PAP'.
        /// </summary>
        void Factorize(CompressedColumnStorage<double> A)
        {
            int n = A.ColumnCount;

            var ap = A.ColumnPointers;
            var ai = A.RowIndices;
            var ax = A.Values;

            int[] parent = S.parent;
            int[] P = S.q;
            int[] Pinv = S.pinv;

            this.D = new double[n];
            this.L = CompressedColumnStorage<double>.Create(n, n, S.cp[n]);

            Array.Copy(S.cp, L.ColumnPointers, n + 1);

            var lp = L.ColumnPointers;
            var li = L.RowIndices;
            var lx = L.Values;

            // Workspace
            var y = new double[n];
            var pattern = new int[n];
            var flag = new int[n];
            var lnz = new int[n];

            double yi, l_ki;
            int i, k, p, kk, p2, len, top;

            for (k = 0; k < n; k++)
            {
                // compute nonzero Pattern of kth row of L, in topological order
                y[k] = 0.0; // Y(0:k) is now all zero
                top = n; // stack for pattern is empty
                flag[k] = k; // mark node k as visited
                lnz[k] = 0; // count of nonzeros in column k of L
                kk = (P != null) ? (P[k]) : (k); // kth original, or permuted, column
                p2 = ap[kk + 1];
                for (p = ap[kk]; p < p2; p++)
                {
                    i = (Pinv != null) ? (Pinv[ai[p]]) : (ai[p]); // get A(i,k)
                    if (i <= k)
                    {
                        y[i] += ax[p]; // scatter A(i,k) into Y (sum duplicates)
                        for (len = 0; flag[i] != k; i = parent[i])
                        {
                            pattern[len++] = i; // L(k,i) is nonzero
                            flag[i] = k; // mark i as visited
                        }
                        while (len > 0)
                        {
                            pattern[--top] = pattern[--len];
                        }
                    }
                }

                // compute numerical values kth row of L (a sparse triangular solve)
                D[k] = y[k]; // get D(k,k) and clear Y(k)
                y[k] = 0.0;
                for (; top < n; top++)
                {
                    i = pattern[top]; // Pattern [top:n-1] is pattern of L(:,k)
                    yi = y[i]; // get and clear Y(i)
                    y[i] = 0.0;
                    p2 = lp[i] + lnz[i];
                    for (p = lp[i]; p < p2; p++)
                    {
                        y[li[p]] -= lx[p] * yi;
                    }
                    l_ki = yi / D[i]; // the nonzero entry L(k,i)
                    D[k] -= l_ki * yi;
                    li[p] = k; // store L(k,i) in column form of L
                    lx[p] = l_ki;
                    lnz[i]++; // increment count of nonzeros in col i
                }

                if (D[k] == 0.0)
                {
                    // failure, D(k,k) is zero
                    throw new Exception("Diagonal element is zero.");
                }
            }
        }
    }
}
