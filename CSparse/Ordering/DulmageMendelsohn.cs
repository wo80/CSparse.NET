namespace CSparse.Ordering
{
    using CSparse.Storage;
    using System;

    /// <summary>
    /// Dulmage-Mendelsohn decomposition.
    /// </summary>
    /// <remarks>
    /// See Chapter 7.4 (Fill-reducing orderings: Dulmage-Mendelsohn decomposition) 
    /// in "Direct Methods for Sparse Linear Systems" by Tim Davis.
    /// </remarks>
    public sealed class DulmageMendelsohn
    {
        private int[] p;   // row permutation (size m)
        private int[] q;   // column permutation (size n)
        private int[] r;   // block k is rows r[k] to r[k+1]-1 in A(p,q) (size nb+1)
        private int[] s;   // block k is cols s[k] to s[k+1]-1 in A(p,q) (size nb+1)

        private int[] rr;  // coarse row decomposition
        private int[] cc;  // coarse column decomposition
        private int nb;    // number of blocks in fine dmperm decomposition

        /// <summary>
        /// Create a new Decomposition instance. 
        /// </summary>
        private DulmageMendelsohn(int m, int n)
        {
            p = new int[m];
            r = new int[m + 6];
            q = new int[n];
            s = new int[n + 6];

            rr = new int[5];  // coarse row decomposition
            cc = new int[5];  // coarse column decomposition
        }

        /// <summary>
        /// Gets the number of blocks in the fine decomposition.
        /// </summary>
        public int Blocks => nb;

        /// <summary>
        /// Gets the number structural rank of the matrix.
        /// </summary>
        public int StructuralRank => rr[3];

        /// <summary>
        /// Gets the number of singletons.
        /// </summary>
        public int Singletons
        {
            get
            {
                int ns = 0;

                for (int k = 0; k < nb; k++)
                {
                    if ((r[k + 1] == r[k] + 1) && (s[k + 1] == s[k] + 1))
                    {
                        ns++;
                    }
                }

                return ns;
            }
        }

        /// <summary>
        /// Gets the row permutation.
        /// </summary>
        public int[] RowPermutation => p;

        /// <summary>
        /// Gets the column permutation.
        /// </summary>
        public int[] ColumnPermutation => q;

        /// <summary>
        /// Gets the block row pointers (block k is rows r[k] to r[k+1]-1 in A(p,q)).
        /// </summary>
        public int[] BlockRowPointers => r;

        /// <summary>
        /// Gets the block column pointers (block k is cols s[k] to s[k+1]-1 in A(p,q)).
        /// </summary>
        public int[] BlockColumnPointers => s;

        /// <summary>
        /// Gets the coarse row decomposition
        /// </summary>
        public int[] CoarseRowDecomposition => rr;

        /// <summary>
        /// Gets the coarse column decomposition.
        /// </summary>
        public int[] CoarseColumnDecomposition => cc;

        /// <summary>
        /// Compute coarse and then fine Dulmage-Mendelsohn decomposition. seed
        /// optionally selects a randomized algorithm.
        /// </summary>
        /// <param name="matrix">column-compressed matrix</param>
        /// <param name="seed">0: natural, -1: reverse, random order otherwise</param>
        /// <returns>Dulmage-Mendelsohn analysis</returns>
        public static DulmageMendelsohn Generate<T>(CompressedColumnStorage<T> matrix, int seed = 0)
             where T : struct, IEquatable<T>, IFormattable
        {
            return Generate(SymbolicColumnStorage.Create(matrix), seed);
        }

        /// <summary>
        /// Compute coarse and then fine Dulmage-Mendelsohn decomposition. seed
        /// optionally selects a randomized algorithm.
        /// </summary>
        /// <param name="A">The matrix represented by <see cref="SymbolicColumnStorage"/>.</param>
        /// <param name="seed">0: natural, -1: reverse, random order otherwise</param>
        /// <returns>Dulmage-Mendelsohn analysis</returns>
        public static DulmageMendelsohn Generate(SymbolicColumnStorage A, int seed = 0)
        {
            int i, j, k, cnz, nc, nb1, nb2;
            int[] Cp, ps, rs;
            bool ok;

            // Maximum matching
            int m = A.RowCount;
            int n = A.ColumnCount;

            var result = new DulmageMendelsohn(m, n); // allocate result
            int[] p = result.p;
            int[] q = result.q;
            int[] r = result.r;
            int[] s = result.s;
            int[] cc = result.cc;
            int[] rr = result.rr;

            int[] jimatch = MaximumMatching.Generate(A, seed); // max transversal

            if (jimatch == null) return null;

            // Coarse decomposition
            for (j = 0; j < n; j++) s[j] = -1; // unmark all cols for bfs
            for (i = 0; i < m; i++) r[i] = -1; // unmark all rows for bfs

            result.BreadthFirstSearch(A, n, r, s, q, jimatch, m, 0, 1); // find C1, R1 from C0*/
            ok = result.BreadthFirstSearch(A, m, s, r, p, jimatch, 0, m, 3); // find R3, C3 from R0*/

            if (!ok) return null;

            result.Unmatched(n, s, q, cc, 0); // unmatched set C0
            result.Matched(n, s, jimatch, m, p, q, cc, rr, 1, 1); // set R1 and C1
            result.Matched(n, s, jimatch, m, p, q, cc, rr, 2, -1); // set R2 and C2
            result.Matched(n, s, jimatch, m, p, q, cc, rr, 3, 3); // set R3 and C3
            result.Unmatched(m, r, p, rr, 3); // unmatched set R0

            // Fine decomposition
            int[] pinv = Permutation.Invert(p); // pinv=p'

            var C = A.Clone();
            A.Permute(pinv, q, C); // C=A(p,q) (it will hold A(R2,C2))

            Cp = C.ColumnPointers;
            nc = cc[3] - cc[2]; // delete cols C0, C1, and C3 from C
            if (cc[2] > 0)
            {
                for (j = cc[2]; j <= cc[3]; j++)
                {
                    Cp[j - cc[2]] = Cp[j];
                }
            }
            C.Reshape(-1, nc);
            if (rr[2] - rr[1] < m) // delete rows R0, R1, and R3 from C
            {
                RowPrune(C, nc, rr);
                cnz = Cp[nc];
                int[] Ci = C.RowIndices;
                if (rr[1] > 0)
                {
                    for (k = 0; k < cnz; k++)
                    {
                        Ci[k] -= rr[1];
                    }
                }
            }
            C.Reshape(nc, -1);
            var scc = StronglyConnectedComponents.Generate(C, nc); // find strongly connected components of C*/

            // Combine coarse and fine decompositions
            ps = scc.p; // C(ps,ps) is the permuted matrix
            rs = scc.r; // kth block is rs[k]..rs[k+1]-1
            nb1 = scc.nb; // # of blocks of A(R2,C2)
            for (k = 0; k < nc; k++) s[k] = q[ps[k] + cc[2]];
            for (k = 0; k < nc; k++) q[k + cc[2]] = s[k];
            for (k = 0; k < nc; k++) r[k] = p[ps[k] + rr[1]];
            for (k = 0; k < nc; k++) p[k + rr[1]] = r[k];
            nb2 = 0; // create the fine block partitions
            r[0] = s[0] = 0;
            if (cc[2] > 0) nb2++; // leading coarse block A (R1, [C0 C1])
            for (k = 0; k < nb1; k++) // coarse block A (R2,C2)
            {
                r[nb2] = rs[k] + rr[1]; // A (R2,C2) splits into nb1 fine blocks
                s[nb2] = rs[k] + cc[2];
                nb2++;
            }
            if (rr[2] < m)
            {
                r[nb2] = rr[2]; // trailing coarse block A ([R3 R0], C3)
                s[nb2] = cc[3];
                nb2++;
            }
            r[nb2] = m;
            s[nb2] = n;
            result.nb = nb2;

            // Remove unused space
            Array.Resize(ref result.r, nb2 + 1);
            Array.Resize(ref result.s, nb2 + 1);

            return result;
        }

        #region Dulmage-Mendelsohn helper

        // breadth-first search for coarse decomposition (C0,C1,R1 or R0,R3,C3)
        private bool BreadthFirstSearch(SymbolicColumnStorage A, int n, int[] wi, int[] wj, int[] queue,
            int[] jimatch, int imatch_offset, int jmatch_offset, int mark)
        {
            // cs_bfs
            int[] Ap, Ai;
            int head = 0, tail = 0, j, i, p, j2;

            for (j = 0; j < n; j++) // place all unmatched nodes in queue
            {
                if (jimatch[imatch_offset + j] >= 0) continue; // skip j if matched
                wj[j] = 0; // j in set C0 (R0 if transpose)
                queue[tail++] = j; // place unmatched col j in queue
            }
            if (tail == 0) return true; // quick return if no unmatched nodes

            // Transpose if requested
            SymbolicColumnStorage C = (mark == 1) ? A.Clone() : A.Transpose();

            if (C == null) return false; // bfs of C=A' to find R3,C3 from R0
            Ap = C.ColumnPointers;
            Ai = C.RowIndices;
            while (head < tail) // while queue is not empty
            {
                j = queue[head++]; // get the head of the queue
                for (p = Ap[j]; p < Ap[j + 1]; p++)
                {
                    i = Ai[p];
                    if (wi[i] >= 0) continue; // skip if i is marked
                    wi[i] = mark; // i in set R1 (C3 if transpose)
                    j2 = jimatch[jmatch_offset + i]; // traverse alternating path to j2
                    if (wj[j2] >= 0) continue; // skip j2 if it is marked
                    wj[j2] = mark; // j2 in set C1 (R3 if transpose)
                    queue[tail++] = j2; // add j2 to queue
                }
            }
            //if (mark != 1) SparseMatrix.spfree(C); // free A' if it was created
            return true;
        }

        // collect matched rows and columns into p and q
        private void Matched(int n, int[] wj, int[] imatch, int imatch_offset, int[] p, int[] q,
            int[] cc, int[] rr, int set, int mark)
        {
            int kc = cc[set];
            int kr = rr[set - 1];
            for (int j = 0; j < n; j++)
            {
                if (wj[j] != mark) continue; // skip if j is not in C set
                p[kr++] = imatch[imatch_offset + j];
                q[kc++] = j;
            }
            cc[set + 1] = kc;
            rr[set] = kr;
        }

        // collect unmatched rows into the permutation vector p
        private void Unmatched(int m, int[] wi, int[] p, int[] rr, int set)
        {
            int i, kr = rr[set];
            for (i = 0; i < m; i++) if (wi[i] == 0) p[kr++] = i;
            rr[set + 1] = kr;
        }

        // TODO: remove unused RowPrune

        // return 1 if row i is in R2
        private static bool RowPrune(int i, int j, int[] rr)
        {
            return (i >= rr[1] && i < rr[2]);
        }

        /// <summary>
        /// Drops entries from a sparse matrix
        /// </summary>
        private static int RowPrune(SymbolicColumnStorage A, int n, int[] rr)
        {
            int i, j, nz = 0;

            var ai = A.RowIndices;
            var ap = A.ColumnPointers;

            for (j = 0; j < n; j++)
            {
                i = ap[j];

                // Record new location of col j.
                ap[j] = nz;

                for (; i < ap[j + 1]; i++)
                {
                    if (ai[i] >= rr[1] && ai[i] < rr[2])
                    {
                        // Keep A(i,j).
                        ai[nz] = ai[i];
                        nz++;
                    }
                }
            }

            // Record new nonzero count.
            ap[n] = nz;

            // Remove extra space.
            Array.Resize(ref A.RowIndices, nz);

            return nz;
        }

        #endregion
    }
}
