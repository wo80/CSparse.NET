namespace CSparse.Factorization
{
    using CSparse.Storage;
    using System;

    /// <summary>
    /// Sparse QR decomposition abstract base class.
    /// </summary>
    public abstract class SparseQR<T> : ISparseFactorization<T>
        where T : struct, IEquatable<T>, IFormattable
    {
        protected readonly int m, n;

        protected SymbolicFactorization S;
        protected CompressedColumnStorage<T> Q, R;
        protected double[] beta;

        /// <summary>
        /// Initializes a new instance of the SparseQR class.
        /// </summary>
        protected SparseQR(int rows, int columns)
        {
            this.m = rows;
            this.n = columns;
        }

        /// <summary>
        /// Gets the number of nonzeros in both Q and R factors together.
        /// </summary>
        public int NonZerosCount
        {
            get { return (Q.NonZerosCount + R.NonZerosCount - m); }
        }

        /// <summary>
        /// Solves a linear system Ax=b.
        /// </summary>
        /// <param name="input">Right hand side b.</param>
        /// <param name="result">Solution vector x.</param>
        public abstract void Solve(T[] input, T[] result);

        /// <summary>
        /// Create a Householder reflection.
        /// </summary>
        protected abstract T CreateHouseholder(T[] x, int offset, ref double beta, int n);

        /// <summary>
        /// Apply the ith Householder vector to x.
        /// </summary>
        protected abstract bool ApplyHouseholder(CompressedColumnStorage<T> V, int i, double beta, T[] x);

        /// <summary>
        /// Sparse QR factorization [V,beta,pinv,R] = qr(A)
        /// </summary>
        /// <param name="A">The matrix to factorize.</param>
        /// <param name="progress">Report progress (range from 0.0 to 1.0).</param>
        protected void Factorize(CompressedColumnStorage<T> A, IProgress<double> progress)
        {
            T zero = Helper.ZeroOf<T>();

            int i, j, p, p1, top, len, col;

            int n = A.ColumnCount;

            var ap = A.ColumnPointers;
            var ai = A.RowIndices;
            var ax = A.Values;

            int[] q = S.q;
            int[] parent = S.parent;
            int[] pinv = S.pinv;
            int m2 = S.m2;

            int vnz = S.lnz;
            int rnz = S.unz;

            int[] leftmost = S.leftmost;

            int[] w = new int[m2 + n]; // get int workspace
            T[] x = new T[m2]; // get double workspace

            int s = m2; // offset into w

            // Allocate result V, R and beta
            var V = this.Q = CompressedColumnStorage<T>.Create(m2, n, vnz);
            var R = this.R = CompressedColumnStorage<T>.Create(m2, n, rnz);
            var b = this.beta = new double[n];

            var rp = R.ColumnPointers;
            var ri = R.RowIndices;
            var rx = R.Values;

            var vp = V.ColumnPointers;
            var vi = V.RowIndices;
            var vx = V.Values;

            for (i = 0; i < m2; i++)
            {
                w[i] = -1; // clear w, to mark nodes
            }

            double current = 0.0;
            double step = n / 100.0;

            rnz = 0; vnz = 0;
            for (int k = 0; k < n; k++) // compute V and R
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

                rp[k] = rnz;      // R(:,k) starts here
                vp[k] = p1 = vnz; // V(:,k) starts here
                w[k] = k;         // add V(k,k) to pattern of V
                vi[vnz++] = k;
                top = n;
                col = q != null ? q[k] : k;
                for (p = ap[col]; p < ap[col + 1]; p++) // find R(:,k) pattern
                {
                    i = leftmost[ai[p]]; // i = min(find(A(i,q)))
                    for (len = 0; w[i] != k; i = parent[i]) // traverse up to k
                    {
                        //len++;
                        w[s + len++] = i;
                        w[i] = k;
                    }
                    while (len > 0)
                    {
                        --top;
                        --len;
                        w[s + top] = w[s + len]; // push path on stack
                    }
                    i = pinv[ai[p]]; // i = permuted row of A(:,col)
                    x[i] = ax[p];    // x (i) = A(:,col)
                    if (i > k && w[i] < k) // pattern of V(:,k) = x (k+1:m)
                    {
                        vi[vnz++] = i; // add i to pattern of V(:,k)
                        w[i] = k;
                    }
                }
                for (p = top; p < n; p++) // for each i in pattern of R(:,k)
                {
                    i = w[s + p]; // R(i,k) is nonzero
                    ApplyHouseholder(V, i, b[i], x); // apply (V(i),Beta(i)) to x
                    ri[rnz] = i; // R(i,k) = x(i)
                    rx[rnz++] = x[i];
                    x[i] = zero;
                    if (parent[i] == k)
                    {
                        //vnz = V.Scatter(i, zero, w, null, k, V, vnz);
                        for (int l = vp[i]; l < vp[i + 1]; l++)
                        {
                            j = vi[l]; // V(j,i) is nonzero
                            if (w[j] < k)
                            {
                                w[j] = k; // j is new entry in column i
                                vi[vnz++] = j; // add j to pattern of V(:,i)
                            }
                        }
                    }
                }
                for (p = p1; p < vnz; p++) // gather V(:,k) = x
                {
                    vx[p] = x[vi[p]];
                    x[vi[p]] = zero;
                }
                ri[rnz] = k; // R(k,k) = norm (x)
                rx[rnz++] = CreateHouseholder(vx, p1, ref b[k], vnz - p1); // [v,beta]=house(x)
            }

            rp[n] = rnz; // finalize R
            vp[n] = vnz; // finalize V
        }

        /// <summary>
        /// Symbolic ordering and analysis for QR.
        /// </summary>
        /// <param name="A">Matrix to factorize.</param>
        /// <param name="p">Permutation.</param>
        /// <param name="natural">Indicates whether to use natural ordering or given permutation.</param>
        protected void SymbolicAnalysis(CompressedColumnStorage<T> A, int[] p, bool natural)
        {
            int m = A.RowCount;
            int n = A.ColumnCount;

            var sym = this.S = new SymbolicFactorization();

            // Fill-reducing ordering
            sym.q = p;

            var C = natural ? SymbolicColumnStorage.Create(A) : Permute(A, null, sym.q);

            // etree of C'*C, where C=A(:,q)
            sym.parent = GraphHelper.EliminationTree(m, n, C.ColumnPointers, C.RowIndices, true);
            int[] post = GraphHelper.TreePostorder(sym.parent, n);
            sym.cp = GraphHelper.ColumnCounts(C, sym.parent, post, true); // col counts chol(C'*C)

            bool ok = C != null && sym.parent != null && sym.cp != null && CountV(C, sym);

            if (ok)
            {
                sym.unz = 0;
                for (int k = 0; k < n; k++)
                {
                    sym.unz += sym.cp[k];
                }
            }
        }

        /// <summary>
        /// Compute nnz(V) = S.lnz, S.pinv, S.leftmost, S.m2 from A and S.parent
        /// </summary>
        private bool CountV(SymbolicColumnStorage A, SymbolicFactorization S)
        {
            int i, k, p, pa;
            int[] pinv, leftmost, parent = S.parent;

            int m = A.RowCount;
            int n = A.ColumnCount;
            int[] ap = A.ColumnPointers;
            int[] ai = A.RowIndices;

            S.pinv = pinv = new int[m + n]; // allocate pinv,
            S.leftmost = leftmost = new int[m]; // and leftmost

            var w = new int[m]; // get workspace
            var head = new int[n];
            var tail = new int[n];
            var nque = new int[n]; // Initialized to 0's

            for (k = 0; k < n; k++) head[k] = -1; // queue k is empty
            for (k = 0; k < n; k++) tail[k] = -1;
            for (i = 0; i < m; i++) leftmost[i] = -1;
            for (k = n - 1; k >= 0; k--)
            {
                for (p = ap[k]; p < ap[k + 1]; p++)
                {
                    leftmost[ai[p]] = k; // leftmost[i] = min(find(A(i,:)))
                }
            }
            for (i = m - 1; i >= 0; i--) // scan rows in reverse order
            {
                pinv[i] = -1; // row i is not yet ordered
                k = leftmost[i];
                if (k == -1) continue; // row i is empty
                if (nque[k]++ == 0) tail[k] = i; // first row in queue k
                w[i] = head[k]; // put i at head of queue k
                head[k] = i;
            }
            S.lnz = 0;
            S.m2 = m;
            for (k = 0; k < n; k++) // find row permutation and nnz(V)
            {
                i = head[k]; // remove row i from queue k
                S.lnz++; // count V(k,k) as nonzero
                if (i < 0) i = S.m2++; // add a fictitious row
                pinv[i] = k; // associate row i with V(:,k)
                if (--nque[k] <= 0) continue; // skip if V(k+1:m,k) is empty
                S.lnz += nque[k]; // nque [k] is nnz (V(k+1:m,k))
                if ((pa = parent[k]) != -1) // move all rows to parent of k
                {
                    if (nque[pa] == 0) tail[pa] = tail[k];
                    w[tail[k]] = head[pa];
                    head[pa] = w[i];
                    nque[pa] += nque[k];
                }
            }
            for (i = 0; i < m; i++)
            {
                if (pinv[i] < 0)
                {
                    pinv[i] = k++;
                }
            }

            return true;
        }

        /// <summary>
        /// Permutes a sparse matrix, C = PAQ.
        /// </summary>
        /// <param name="A">m-by-n, column-compressed matrix</param>
        /// <param name="pinv">a permutation vector of length m</param>
        /// <param name="q">a permutation vector of length n</param>
        /// <returns>C = PAQ, null on error</returns>
        private SymbolicColumnStorage Permute(CompressedColumnStorage<T> A, int[] pinv, int[] q)
        {
            int t, j, k, nz = 0;

            int m = A.RowCount;
            int n = A.ColumnCount;

            var ap = A.ColumnPointers;
            var ai = A.RowIndices;

            var result = SymbolicColumnStorage.Create(A);

            var cp = result.ColumnPointers;
            var ci = result.RowIndices;

            for (k = 0; k < n; k++)
            {
                // Column k of C is column q[k] of A
                cp[k] = nz;
                j = q != null ? (q[k]) : k;
                for (t = ap[j]; t < ap[j + 1]; t++)
                {
                    // Row i of A is row pinv[i] of C
                    ci[nz++] = pinv != null ? (pinv[ai[t]]) : ai[t];
                }
            }

            // Finalize the last column of C
            cp[n] = nz;

            return result;
        }
    }
}
