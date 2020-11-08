namespace CSparse
{
    using CSparse.Storage;
    using System;
    using System.Collections.Generic;

    /// <summary>
    /// Converter for different types of storages.
    /// </summary>
    public static class Converter
    {
        /// <summary>
        /// Convert a coordinate storage to compressed sparse column (CSC) format.
        /// </summary>
        /// <param name="storage">Coordinate storage.</param>
        /// <param name="cleanup">Remove and sum duplicate entries.</param>
        /// <returns>Compressed sparse column storage.</returns>
        [Obsolete("Will be removed in future versions. Use SparseMatrix.OfIndexed(...) instead.")]
        public static CompressedColumnStorage<T> ToCompressedColumnStorage<T>(CoordinateStorage<T> storage,
            bool cleanup = true) where T : struct, IEquatable<T>, IFormattable
        {
            return ToCompressedColumnStorage_(storage, cleanup);
        }

        /// <summary>
        /// Convert a coordinate storage to compressed sparse column (CSC) format.
        /// </summary>
        /// <param name="storage">Coordinate storage.</param>
        /// <param name="cleanup">Remove and sum duplicate entries.</param>
        /// <param name="inplace">Do the conversion in place (re-using the coordinate storage arrays).</param>
        /// <returns>Compressed sparse column storage.</returns>
        internal static CompressedColumnStorage<T> ToCompressedColumnStorage_<T>(CoordinateStorage<T> storage,
                bool cleanup = true, bool inplace = false) where T : struct, IEquatable<T>, IFormattable
        {
            int nrows = storage.RowCount;
            int ncols = storage.ColumnCount;

            int nz = storage.NonZerosCount;

            var result = CompressedColumnStorage<T>.Create(nrows, ncols);

            var ap = result.ColumnPointers = new int[ncols + 1];

            if (nz == 0)
            {
                return result;
            }

            var values = storage.Values;
            var rowind = storage.RowIndices;
            var colind = storage.ColumnIndices;

            if (inplace)
            {
                ConvertInPlace(ncols, nz, values, rowind, colind, ap);

                result.RowIndices = rowind;
                result.Values = values;

                // Make sure the data can't be accessed through the coordinate storage.
                storage.Invalidate();
            }
            else
            {
                var columnCounts = new int[ncols];

                for (int k = 0; k < nz; k++)
                {
                    // Count columns
                    columnCounts[colind[k]]++;
                }

                // Get column pointers
                int valueCount = Helper.CumulativeSum(ap, columnCounts, ncols);

                var ai = new int[valueCount];
                var ax = new T[valueCount];

                for (int k = 0; k < nz; k++)
                {
                    int p = columnCounts[colind[k]]++;
                    ai[p] = rowind[k];
                    ax[p] = values[k];
                }

                result.RowIndices = ai;
                result.Values = ax;
            }

            Helper.SortIndices(result);

            if (cleanup)
            {
                result.Cleanup();
            }

            return result;
        }

        /// <summary>
        /// Converts a matrix stored in coordinate format into the compressed sparse
        /// column (CSC) format. The conversion is done in place.
        /// </summary>
        /// <param name="columns">Number of columns.</param>
        /// <param name="nz">Number of nonzero elements.</param>
        /// <param name="values">Nonzero elements (size nz).</param>
        /// <param name="colind">Column indices (size nz).</param>
        /// <param name="rowind">Row indices (size nz).</param>
        /// <param name="work">Work array (size n+1).</param>
        /// <remarks>
        /// On return, the coordinate storage input arrays contain the compressed sparse
        /// column data structure for the resulting matrix. The <paramref name="work"/>
        /// array contains a copy of the column pointer.
        /// 
        /// The entries of the output matrix are not sorted (the row indices in each
        /// column are not in increasing order).
        /// </remarks>
        private static void ConvertInPlace<T>(int columns, int nz, T[] values, int[] rowind, int[] colind, int[] work)
            where T : struct, IEquatable<T>, IFormattable
        {
            int i, j;
            int pos, inext, jnext;
            T t = default;
            T tnext = default;

            // Count column entries.
            for (i = 0; i < nz; i++)
            {
                work[colind[i] + 1]++;
            }

            // Cumulative sum.
            work[0] = 0;
            for (i = 1; i < columns + 1; i++)
            {
                work[i] = work[i - 1] + work[i];
            }

            int init = 0;
            int k = 0;

            // Start chasing process.
            while (init < nz)
            {
                t = values[init];
                i = rowind[init];
                j = colind[init];

                colind[init] = -1;

                while (k < nz)
                {
                    k++;

                    // Current column number is j. Determine where to go.
                    pos = work[j];

                    // Save the chased element.
                    tnext = values[pos];
                    inext = rowind[pos];
                    jnext = colind[pos];

                    // Then occupy its location.
                    values[pos] = t;
                    rowind[pos] = i;

                    // Update pointer information for next element to come in column j.
                    work[j] = pos + 1;

                    if (colind[pos] < 0) break;

                    // Determine next element to be chased.
                    t = tnext;
                    i = inext;
                    j = jnext;

                    colind[pos] = -1;
                }

                while (++init < nz && colind[init] < 0) ;

                // Restart chasing.
            }

            // Fix column pointers.
            Buffer.BlockCopy(work, 0, work, Constants.SizeOfInt, columns * Constants.SizeOfInt);

            work[0] = 0;
        }

        /// <summary>
        /// Convert a 2D array to coordinate storage.
        /// </summary>
        /// <param name="array">2D array storage.</param>
        /// <returns>Coordinate storage.</returns>
        [Obsolete("Will be removed in future versions. Use SparseMatrix.OfArray(...) instead.")]
        public static CoordinateStorage<T> FromDenseArray<T>(T[,] array)
            where T : struct, IEquatable<T>, IFormattable
        {
            return FromDenseArray_(array);
        }

        internal static CoordinateStorage<T> FromDenseArray_<T>(T[,] array)
            where T : struct, IEquatable<T>, IFormattable
        {
            int rowCount = array.GetLength(0);
            int columnCount = array.GetLength(1);

            var storage = new CoordinateStorage<T>(rowCount, columnCount, Math.Max(rowCount, columnCount));

            for (int i = 0; i < rowCount; i++)
            {
                for (int j = 0; j < columnCount; j++)
                {
                    storage.At(i, j, array[i, j]);
                }
            }

            return storage;
        }

        /// <summary>
        /// Convert a jagged array to compressed sparse column (CSC) format.
        /// </summary>
        /// <param name="array">Jagged array storage.</param>
        /// <returns>Compressed sparse column storage.</returns>
        [Obsolete("Will be removed in future versions. Use SparseMatrix.OfJaggedArray(...) instead.")]
        public static CompressedColumnStorage<T> ToCompressedColumnStorage<T>(T[][] array)
            where T : struct, IEquatable<T>, IFormattable
        {
            int nrows = array.Length;
            int ncols = array[0].Length;

            var storage = new CoordinateStorage<T>(nrows, ncols, nrows);

            for (int i = 0; i < nrows; i++)
            {
                for (int j = 0; j < ncols; j++)
                {
                    storage.At(i, j, array[i][j]);
                }
            }

            return ToCompressedColumnStorage_<T>(storage, false);
        }


        /// <summary>
        /// Convert a column major array to coordinate storage.
        /// </summary>
        /// <param name="array">Column major array storage.</param>
        /// <param name="rowCount">Number of rows.</param>
        /// <param name="columnCount">Number of columns.</param>
        /// <returns>Coordinate storage.</returns>
        [Obsolete("Will be removed in future versions. Use SparseMatrix.OfColumnMajor(...) instead.")]
        public static CoordinateStorage<T> FromColumnMajorArray<T>(T[] array, int rowCount, int columnCount)
            where T : struct, IEquatable<T>, IFormattable
        {
            return FromColumnMajorArray_(array, rowCount, columnCount);
        }

        internal static CoordinateStorage<T> FromColumnMajorArray_<T>(T[] array, int rowCount, int columnCount)
            where T : struct, IEquatable<T>, IFormattable
        {
            var storage = new CoordinateStorage<T>(rowCount, columnCount, Math.Max(rowCount, columnCount));

            for (int i = 0; i < rowCount; i++)
            {
                for (int j = 0; j < columnCount; j++)
                {
                    storage.At(i, j, array[i + j * rowCount]);
                }
            }

            return storage;
        }

        /// <summary>
        /// Convert a 2D jagged array to coordinate storage.
        /// </summary>
        /// <param name="array">jagged array storage.</param>
        /// <returns>Coordinate storage.</returns>
        /// <remarks>All rows of the array are assumed to be equal in length</remarks>
        [Obsolete("Will be removed in future versions. Use SparseMatrix.OfColumnMajor(...) instead.")]
        public static CoordinateStorage<T> FromJaggedArray<T>(T[][] array)
            where T : struct, IEquatable<T>, IFormattable
        {
            return FromJaggedArray_(array);
        }

        internal static CoordinateStorage<T> FromJaggedArray_<T>(T[][] array)
            where T : struct, IEquatable<T>, IFormattable
        {
            int rowCount = array.Length;
            int columnCount = array[0].Length;

            var storage = new CoordinateStorage<T>(rowCount, columnCount, rowCount);

            for (int i = 0; i < rowCount; i++)
            {
                for (int j = 0; j < columnCount; j++)
                {
                    storage.At(i, j, array[i][j]);
                }
            }

            return storage;
        }

        /// <summary>
        /// Convert a row major array to coordinate storage.
        /// </summary>
        /// <param name="array">Row major array storage.</param>
        /// <param name="rowCount">Number of rows.</param>
        /// <param name="columnCount">Number of columns.</param>
        /// <returns>Coordinate storage.</returns>
        [Obsolete("Will be removed in future versions. Use SparseMatrix.OfRowMajor(...) instead.")]
        public static CoordinateStorage<T> FromRowMajorArray<T>(T[] array, int rowCount, int columnCount)
            where T : struct, IEquatable<T>, IFormattable
        {
            return FromRowMajorArray_(array, rowCount, columnCount);
        }

        internal static CoordinateStorage<T> FromRowMajorArray_<T>(T[] array, int rowCount, int columnCount)
            where T : struct, IEquatable<T>, IFormattable
        {
            var storage = new CoordinateStorage<T>(rowCount, columnCount, Math.Max(rowCount, columnCount));

            for (int i = 0; i < rowCount; i++)
            {
                for (int j = 0; j < columnCount; j++)
                {
                    storage.At(i, j, array[i * columnCount + j]);
                }
            }

            return storage;
        }

        /// <summary>
        /// Convert a row major array to coordinate storage.
        /// </summary>
        /// <param name="enumerable">Enumerates the entries of a matrix.</param>
        /// <param name="rowCount">Number of rows.</param>
        /// <param name="columnCount">Number of columns.</param>
        /// <returns>Coordinate storage.</returns>
        [Obsolete("Will be removed in future versions. Use SparseMatrix.OfIndexed(...) instead.")]
        public static CoordinateStorage<T> FromEnumerable<T>(IEnumerable<Tuple<int, int, T>> enumerable, int rowCount, int columnCount)
            where T : struct, IEquatable<T>, IFormattable
        {
            return FromEnumerable_(enumerable, rowCount, columnCount);
        }

        internal static CoordinateStorage<T> FromEnumerable_<T>(IEnumerable<Tuple<int, int, T>> enumerable, int rowCount, int columnCount)
            where T : struct, IEquatable<T>, IFormattable
        {
            var storage = new CoordinateStorage<T>(rowCount, columnCount, Math.Max(rowCount, columnCount));

            foreach (var item in enumerable)
            {
                storage.At(item.Item1, item.Item2, item.Item3);
            }

            return storage;
        }
    }
}
