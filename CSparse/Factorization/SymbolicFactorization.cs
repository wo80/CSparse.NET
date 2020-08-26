namespace CSparse.Factorization
{
    /// <summary>
    /// Symbolic Cholesky, LU, or QR factorization storage.
    /// </summary>
    public class SymbolicFactorization
    {
        /// <summary>
        /// inverse row perm. for QR, fill red. perm for Chol
        /// </summary>
        public int[] pinv;

        /// <summary>
        /// fill-reducing column permutation for LU and QR
        /// </summary>
        public int[] q;

        /// <summary>
        /// elimination tree for Cholesky and QR
        /// </summary>
        public int[] parent;

        /// <summary>
        /// column pointers for Cholesky, row counts for QR
        /// </summary>
        public int[] cp;

        /// <summary>
        /// leftmost[i] = min(find(A(i,:))), for QR
        /// </summary>
        public int[] leftmost;

        /// <summary>
        /// # of rows for QR, after adding fictitious rows
        /// </summary>
        public int m2;

        /// <summary>
        /// # entries in L for LU or Cholesky; in V for QR
        /// </summary>
        public int lnz;

        /// <summary>
        /// # entries in U for LU; in R for QR
        /// </summary>
        public int unz;
    }
}
