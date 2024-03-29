namespace CSparse.Complex
{
    using CSparse.Storage;
    using System;
    using System.Numerics;

    /// <summary>
    /// Helper methods for solving triangular systems.
    /// </summary>
    public static class SolverHelper
    {
        /// <summary>
        /// Solve a lower triangular system by forward elimination, Lx=b.
        /// </summary>
        /// <param name="L"></param>
        /// <param name="x"></param>
        /// <returns></returns>
        public static void SolveLower(CompressedColumnStorage<Complex> L, Complex[] x) => SolveLower(L, x.AsSpan());

        /// <summary>
        /// Solve a lower triangular system by forward elimination, Lx=b.
        /// </summary>
        /// <param name="L"></param>
        /// <param name="x"></param>
        /// <returns></returns>
        public static void SolveLower(CompressedColumnStorage<Complex> L, Span<Complex> x)
        {
            int p, j, k, n = L.ColumnCount;

            var lp = L.ColumnPointers;
            var li = L.RowIndices;
            var lx = L.Values;

            for (j = 0; j < n; j++)
            {
                x[j] /= lx[lp[j]];

                k = lp[j + 1];

                for (p = lp[j] + 1; p < k; p++)
                {
                    x[li[p]] -= lx[p] * x[j];
                }
            }
        }

        /// <summary>
        /// Solve L'x=b where x and b are dense.
        /// </summary>
        /// <param name="L"></param>
        /// <param name="x"></param>
        /// <returns></returns>
        public static void SolveLowerTranspose(CompressedColumnStorage<Complex> L, Complex[] x) => SolveLowerTranspose(L, x.AsSpan());

        /// <summary>
        /// Solve L'x=b where x and b are dense.
        /// </summary>
        /// <param name="L"></param>
        /// <param name="x"></param>
        /// <returns></returns>
        public static void SolveLowerTranspose(CompressedColumnStorage<Complex> L, Span<Complex> x)
        {
            int p, j, k, n = L.ColumnCount;

            var lp = L.ColumnPointers;
            var li = L.RowIndices;
            var lx = L.Values;

            for (j = n - 1; j >= 0; j--)
            {
                k = lp[j + 1];

                for (p = lp[j] + 1; p < k; p++)
                {
                    x[j] -= Complex.Conjugate(lx[p]) * x[li[p]];
                }

                x[j] /= Complex.Conjugate(lx[lp[j]]);
            }
        }

        /// <summary>
        /// Solve an upper triangular system by backward elimination, Ux=b.
        /// </summary>
        /// <param name="U"></param>
        /// <param name="x"></param>
        /// <returns></returns>
        public static void SolveUpper(CompressedColumnStorage<Complex> U, Complex[] x)=>SolveUpper(U, x.AsSpan());

        /// <summary>
        /// Solve an upper triangular system by backward elimination, Ux=b.
        /// </summary>
        /// <param name="U"></param>
        /// <param name="x"></param>
        /// <returns></returns>
        public static void SolveUpper(CompressedColumnStorage<Complex> U, Span<Complex> x)
        {
            int p, j, k, n = U.ColumnCount;

            var up = U.ColumnPointers;
            var ui = U.RowIndices;
            var ux = U.Values;

            for (j = n - 1; j >= 0; j--)
            {
                x[j] /= ux[up[j + 1] - 1];

                k = up[j + 1] - 1;

                for (p = up[j]; p < k; p++)
                {
                    x[ui[p]] -= ux[p] * x[j];
                }
            }
        }

        /// <summary>
        /// Solve U'x=b where x and b are dense.
        /// </summary>
        /// <param name="U"></param>
        /// <param name="x"></param>
        /// <returns></returns>
        public static void SolveUpperTranspose(CompressedColumnStorage<Complex> U, Complex[] x) => SolveUpperTranspose(U, x.AsSpan());

        /// <summary>
        /// Solve U'x=b where x and b are dense.
        /// </summary>
        /// <param name="U"></param>
        /// <param name="x"></param>
        /// <returns></returns>
        public static void SolveUpperTranspose(CompressedColumnStorage<Complex> U, Span<Complex> x)
        {
            int p, j, k, n = U.ColumnCount;

            var up = U.ColumnPointers;
            var ui = U.RowIndices;
            var ux = U.Values;

            for (j = 0; j < n; j++)
            {
                k = up[j + 1] - 1;

                for (p = up[j]; p < k; p++)
                {
                    x[j] -= Complex.Conjugate(ux[p]) * x[ui[p]];
                }

                x[j] /= Complex.Conjugate(ux[up[j + 1] - 1]);
            }
        }
    }
}
