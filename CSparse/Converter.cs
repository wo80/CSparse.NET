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

        internal static CompressedColumnStorage<T> ToCompressedColumnStorage_<T>(CoordinateStorage<T> storage,
                bool cleanup = true) where T : struct, IEquatable<T>, IFormattable
        {
            int nrows = storage.RowCount;
            int ncols = storage.ColumnCount;

            var values = storage.Values;
            var rowind = storage.RowIndices;
            var colind = storage.ColumnIndices;

            int p, k, nz = storage.NonZerosCount;

            var columnPointers = new int[ncols + 1];
            var columnCounts = new int[ncols];

            for (k = 0; k < nz; k++)
            {
                // Count columns
                columnCounts[colind[k]]++;
            }

            // Get row pointers
            int valueCount = Helper.CumulativeSum(columnPointers, columnCounts, ncols);

            var result = CompressedColumnStorage<T>.Create(nrows, ncols);

            var rowIndices = new int[valueCount];
            var storageValues = new T[valueCount];

            for (k = 0; k < nz; k++)
            {
                p = columnCounts[colind[k]]++;
                rowIndices[p] = rowind[k];
                storageValues[p] = values[k];
            }

            result.RowIndices = rowIndices;
            result.ColumnPointers = columnPointers;
            result.Values = storageValues;

            Helper.SortIndices(result);

            if (cleanup)
            {
                result.Cleanup();
            }

            return result;
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
