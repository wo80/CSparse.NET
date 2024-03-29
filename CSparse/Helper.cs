namespace CSparse
{
    using CSparse.Storage;
    using System;

    /// <summary>
    /// Public helper methods.
    /// </summary>
    public static class Helper
    {
        /// <summary>
        /// Cumulative sum of given array.
        /// </summary>
        /// <param name="sum">Output: cumulative sum of counts</param>
        /// <param name="counts">input array, overwritten with sum</param>
        /// <param name="size">length of counts</param>
        /// <returns>sum[size] (non-zeros)</returns>
        public static int CumulativeSum(int[] sum, int[] counts, int size)
        {
            int i, nz = 0;

            for (i = 0; i < size; i++)
            {
                sum[i] = nz;
                nz += counts[i];
                counts[i] = sum[i]; // also copy p[0..n-1] back into c[0..n-1]
            }

            sum[size] = nz;

            return nz;
        }

        /// <summary>
        /// Validate the structure of the <see cref="CompressedColumnStorage{T}"/>.
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="storage">The storage to validate.</param>
        /// <param name="strict">If true, row indices have to be ordered and no duplicate indices are allowed (default = <c>false</c>).</param>
        /// <returns>Returns true if the structure of the storage is valid.</returns>
        public static bool ValidateStorage<T>(CompressedColumnStorage<T> storage, bool strict = false)
            where T : struct, IEquatable<T>, IFormattable
        {
            int rows = storage.RowCount;
            int columns = storage.ColumnCount;

            var ap = storage.ColumnPointers;
            var ai = storage.RowIndices;

            for (int i = 0; i < columns; i++)
            {
                int j = ap[i];
                int end = ap[i + 1];

                // Check if column pointers are in ascending order.
                if (j > end)
                {
                    return false;
                }

                for (; j < end; j++)
                {
                    // Check if row indices are within bounds.
                    if (ai[j] < 0 || ai[j] >= rows)
                    {
                        return false;
                    }

                    // Check if row indices are in order.
                    if (strict && ai[j] >= ai[j + 1])
                    {
                        return false;
                    }
                }
            }

            return true;
        }

        /// <summary>
        /// Trim row indices and values array of the storage to the exact size (non-zeros count).
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="storage">The sparse matrix.</param>
        public static void TrimStorage<T>(CompressedColumnStorage<T> storage)
            where T : struct, IEquatable<T>, IFormattable
        {
            storage.Resize(0);
        }

        /// <summary>
        /// Sort the row indices of the storage.
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="storage">The sparse matrix.</param>
        public static void SortIndices<T>(CompressedColumnStorage<T> storage)
            where T : struct, IEquatable<T>, IFormattable
        {
            int from, size, to, p, q, idx;
            T val;

            int columns = storage.ColumnCount;

            var ap = storage.ColumnPointers;
            var ai = storage.RowIndices;
            var ax = storage.Values;

            for (int i = 0; i < columns; i++)
            {
                from = ap[i];
                to = ap[i + 1] - 1;

                size = to - from + 1;

                if (size > 16)
                {
                    // Quicksort
                    Array.Sort(ai, ax, from, size);
                }
                else
                {
                    // Insertion sort
                    for (p = from + 1; p <= to; p++)
                    {
                        idx = ai[p];
                        val = ax[p];
                        q = p - 1;
                        while (q >= from && ai[q] > idx)
                        {
                            ai[q + 1] = ai[q];
                            ax[q + 1] = ax[q];
                            q--;
                        }
                        ai[q + 1] = idx;
                        ax[q + 1] = val;
                    }
                }
            }
        }

        #region Generic code helper

        /// <summary>
        /// Sets the value of <c>1.0</c> for type T.
        /// </summary>
        /// <typeparam name="T">The type to return the value of 1.0 of.</typeparam>
        /// <returns>The value of <c>1.0</c> for type T.</returns>
        public static T OneOf<T>() where T : struct
        {
            if (typeof(T) == typeof(double))
            {
                return (T)(object)1.0d;
            }

            if (typeof(T) == typeof(System.Numerics.Complex))
            {
                return (T)(object)System.Numerics.Complex.One;
            }

            throw new NotSupportedException();
        }

        /// <summary>
        /// Sets the value of <c>0.0</c> for type T.
        /// </summary>
        /// <typeparam name="T">The type to return the value of 0.0 of.</typeparam>
        /// <returns>The value of <c>0.0</c> for type T.</returns>
        public static T ZeroOf<T>() where T : struct
        {
            return default;
        }

        #endregion
    }
}
