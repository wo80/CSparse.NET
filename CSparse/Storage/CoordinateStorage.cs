// -----------------------------------------------------------------------
// <copyright file="CoordinateStorage.cs">
// Copyright (c) 2012-2016, Christian Woltering
// </copyright>
// -----------------------------------------------------------------------

namespace CSparse.Storage
{
    using System;

    /// <summary>
    /// Coordinate storage sparse matrix format.
    /// </summary>
    public class CoordinateStorage<T>
        where T : struct, IEquatable<T>, IFormattable
    {
        private static readonly T Zero = Helper.ZeroOf<T>();

        private readonly int nrows;
        private readonly int ncols;

        private int nz; // Number of entries in triplet matrix
        private int nzmax; // Maximum number of entries

        private int[] rowind; // Row indices (size nzmax)
        private int[] colind; // Column indices (size nzmax)
        private T[] values; // Numerical values (size nzmax)

        /// <summary>
        /// Row indices (size = NonZerosCount)
        /// </summary>
        public int[] RowIndices
        {
            get { return rowind; }
        }

        /// <summary>
        /// Column indices (size = NonZerosCount)
        /// </summary>
        public int[] ColumnIndices
        {
            get { return colind; }
        }

        /// <summary>
        /// Numerical values (size = NonZerosCount)
        /// </summary>
        public T[] Values
        {
            get { return values; }
        }

        /// <summary>
        /// Gets the number of rows.
        /// </summary>
        public int RowCount
        {
            get { return nrows; }
        }

        /// <summary>
        /// Gets the number of columns.
        /// </summary>
        public int ColumnCount
        {
            get { return ncols; }
        }

        /// <summary>
        /// Gets the number of non-zero entries.
        /// </summary>
        public int NonZerosCount
        {
            get { return nz; }
        }

        /// <summary>
        /// Initializes a new instance of the CoordinateStorage class.
        /// </summary>
        public CoordinateStorage(int rowCount, int columnCount, int nzmax, bool alloc = true)
        {
            this.nrows = rowCount;
            this.ncols = columnCount;

            if (alloc)
            {
                this.nzmax = (nzmax = Math.Max(nzmax, 1));
                this.nz = 0;

                this.rowind = new int[nzmax];
                this.colind = new int[nzmax];
                this.values = new T[nzmax];
            }
        }

        /*
        public void SetStorage(int[] rowIndices, int[] columnIndices, T[] values)
        {
            this.nzmax = values.Length;
            this.nz = values.Length;

            this.rowind = rowIndices;
            this.colind = columnIndices;
            this.values = values;
        }
        //*/

        /// <summary>
        /// Adds an entry. Memory and dimension of the matrix are increased if necessary.
        /// </summary>
        /// <param name="i">Row index of new entry</param>
        /// <param name="j">Column index of new entry</param>
        /// <param name="value">Numerical value of new entry</param>
        /// <returns>True if successful, false otherwise</returns>
        public void At(int i, int j, T value)
        {
            if (i < 0 || j < 0)
            {
                return;
            }

            if (value.Equals(Zero))
            {
                return;
            }

            if (nz >= nzmax)
            {
                Resize(2 * nzmax);
            }

            if (i < 0 || i >= nrows)
            {
                throw new ArgumentOutOfRangeException("i");
            }

            if (j < 0 || j >= ncols)
            {
                throw new ArgumentOutOfRangeException("j");
            }

            rowind[nz] = i;
            colind[nz] = j;
            values[nz] = value;

            nz += 1;
        }

        /// <summary>
        /// Resize the storage arrays of the sparse matrix.
        /// </summary>
        /// <param name="size">The new size of Values and ColumnIndices arrays.</param>
        /// <remarks>
        /// Use size = 0 to automatically resize to non-zeros count.
        /// </remarks>
        protected bool Resize(int size)
        {
            if (size <= 0)
            {
                size = this.nz;
            }

            // TODO: check available memory
            //       and throw OutOfMemoryException before trying to resize.

            Array.Resize(ref this.rowind, size);
            Array.Resize(ref this.colind, size);
            Array.Resize(ref this.values, size);

            this.nzmax = size;

            return true;
        }
    }
}
