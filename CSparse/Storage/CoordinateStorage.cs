namespace CSparse.Storage
{
    using CSparse.Properties;
    using System;

    /// <summary>
    /// Coordinate storage sparse matrix format.
    /// </summary>
    /// <typeparam name="T">Supported data types are <c>double</c> and <see cref="System.Numerics.Complex"/>.</typeparam>
    public class CoordinateStorage<T>
        where T : struct, IEquatable<T>, IFormattable
    {
        private static readonly T Zero = default; // default = zero

        private readonly int nrows;
        private readonly int ncols;

        private int nz; // Number of entries in triplet matrix
        private int nzmax; // Maximum number of entries

        private int[] rowind; // Row indices (size nzmax)
        private int[] colind; // Column indices (size nzmax)
        private T[] values; // Numerical values (size nzmax)

        /// <summary>
        /// Row indices (size = nzmax, may be larger than non-zeros count).
        /// </summary>
        public int[] RowIndices => rowind;

        /// <summary>
        /// Column indices (size = nzmax, may be larger than non-zeros count).
        /// </summary>
        public int[] ColumnIndices => colind;

        /// <summary>
        /// Numerical values (size = nzmax, may be larger than non-zeros count).
        /// </summary>
        public T[] Values => values;

        /// <summary>
        /// Gets the number of rows.
        /// </summary>
        public int RowCount => nrows;

        /// <summary>
        /// Gets the number of columns.
        /// </summary>
        public int ColumnCount => ncols;

        /// <summary>
        /// Gets the number of non-zero entries.
        /// </summary>
        public int NonZerosCount => nz;

        /// <summary>
        /// Initializes a new instance of the <see cref="CoordinateStorage{T}"/> class.
        /// </summary>
        /// <param name="rowCount">The number of rows.</param>
        /// <param name="columnCount">The number of columns.</param>
        /// <param name="nzmax">The number of non-zeros to allocate (will be expanded dynamically).</param>
        public CoordinateStorage(int rowCount, int columnCount, int nzmax)
            : this(rowCount, columnCount, nzmax, nzmax >= 0)
        {
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="CoordinateStorage{T}"/> class.
        /// </summary>
        /// <param name="rowCount">The number of rows.</param>
        /// <param name="columnCount">The number of columns.</param>
        /// <param name="rowind">The row indices array.</param>
        /// <param name="colind">The column indices array.</param>
        /// <param name="values">The values array.</param>
        /// <remarks>
        /// The storage arrays are considered to be empty (non-zeros count is <c>0</c>). If the arrays are
        /// already filled, the constructor taking the explicit <c>nonZerosCount</c> has to be used.
        /// </remarks>
        public CoordinateStorage(int rowCount, int columnCount, int[] rowind, int[] colind, T[] values)
            : this(rowCount, columnCount, 0, rowind, colind, values)
        {
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="CoordinateStorage{T}"/> class.
        /// </summary>
        /// <param name="rowCount">The number of rows.</param>
        /// <param name="columnCount">The number of columns.</param>
        /// <param name="nonZerosCount">The number of non-zeros in the given arrays.</param>
        /// <param name="rowind">The row indices array.</param>
        /// <param name="colind">The column indices array.</param>
        /// <param name="values">The values array.</param>
        public CoordinateStorage(int rowCount, int columnCount, int nonZerosCount, int[] rowind, int[] colind, T[] values)
            : this(rowCount, columnCount, 0, false)
        {
            if (rowind == null)
            {
                throw new ArgumentNullException(nameof(rowind));
            }

            if (colind == null)
            {
                throw new ArgumentNullException(nameof(colind));
            }

            if (values == null)
            {
                throw new ArgumentNullException(nameof(values));
            }

            this.rowind = rowind;
            this.colind = colind;
            this.values = values;

            // Allow arrays to have different sizes.
            nzmax = Math.Min(values.Length, Math.Min(rowind.Length, colind.Length));

            if (nzmax < nonZerosCount)
            {
                throw new ArgumentException("Invalid array storage (all arrays must be at least of size nonZerosCount).");
            }

            nz = nonZerosCount;
        }

        private CoordinateStorage(int rowCount, int columnCount, int nzmax, bool alloc)
        {
            if (rowCount < 0)
            {
                throw new ArgumentOutOfRangeException(nameof(rowCount), Resources.MatrixDimensionNonNegative);
            }

            if (columnCount < 0)
            {
                throw new ArgumentOutOfRangeException(nameof(columnCount), Resources.MatrixDimensionNonNegative);
            }

            if (nzmax < 0)
            {
                throw new ArgumentOutOfRangeException(nameof(nzmax), Resources.ValueNonNegative);
            }

            this.nrows = rowCount;
            this.ncols = columnCount;

            if (alloc)
            {
                this.nzmax = nzmax = Math.Max(nzmax, 1);
                this.nz = 0;

                this.rowind = new int[nzmax];
                this.colind = new int[nzmax];
                this.values = new T[nzmax];
            }
        }

        /// <summary>
        /// Adds an entry to the storage.
        /// </summary>
        /// <param name="i">Row index of new entry</param>
        /// <param name="j">Column index of new entry</param>
        /// <param name="value">Numerical value of new entry</param>
        /// <remarks>
        /// Duplicate entries will be added up, i.e. calling
        /// <code>
        /// storage.At(0, 0, 1.0);
        /// storage.At(0, 0, 2.0);
        /// </code>
        /// will result in an entry with value 3.0 at index (0, 0) of the
        /// resulting matrix.
        /// 
        /// Memory will be increased as necessary.
        /// </remarks>
        public void At(int i, int j, T value)
        {
            if (i < 0 || i >= nrows)
            {
                throw new ArgumentOutOfRangeException(nameof(i));
            }

            if (j < 0 || j >= ncols)
            {
                throw new ArgumentOutOfRangeException(nameof(j));
            }

            if (value.Equals(Zero))
            {
                return;
            }

            if (nz >= nzmax)
            {
                Resize(2 * nzmax);
            }

            rowind[nz] = i;
            colind[nz] = j;
            values[nz] = value;

            nz++;
        }

        /// <summary>
        /// Filter storage values.
        /// </summary>
        /// <param name="func">Filter function returning true if value should be kept,
        /// false if value should be discarded.</param>
        /// <returns>New number of non-zeros.</returns>
        /// <remarks>
        /// Filter function arguments:
        /// 
        /// 1 = Row index i
        /// 2 = Column index j
        /// 3 = Value of entry (i,j)
        /// 
        /// Element a_{i,j} is dropped, if func(i, j, aij) returns false.
        /// </remarks>
        public int Keep(Func<int, int, T, bool> func)
        {
            int k = 0;

            for (int i = 0; i < nz; i++)
            {
                int ai = rowind[i];
                int aj = colind[i];
                var ax = values[i];

                if (func(ai, aj, ax))
                {
                    // Keep A(i,j).
                    rowind[k] = ai;
                    colind[k] = aj;
                    values[k] = ax;
                    k++;
                }
            }

            return nz = k;
        }

        /// <summary>
        /// Remove all values from the storage (without freeing the memory).
        /// </summary>
        public void Clear()
        {
            Array.Clear(rowind, 0, nzmax);
            Array.Clear(colind, 0, nzmax);
            Array.Clear(values, 0, nzmax);

            nz = 0;
        }

        /// <summary>
        /// Returns the transposed coordinate storage.
        /// </summary>
        /// <param name="alloc">If true, clone storage arrays, otherwise just swap the references and re-use the arrays.</param>
        /// <returns>The transposed storage.</returns>
        public CoordinateStorage<T> Transpose(bool alloc = false)
        {
            var result = new CoordinateStorage<T>(ncols, nrows, nzmax, alloc);

            result.nz = nz;
            result.nzmax = nzmax;

            // Transposing is just a matter of switching row and column indices.
            if (alloc)
            {
                Array.Copy(rowind, result.colind, nz);
                Array.Copy(colind, result.rowind, nz);
                Array.Copy(values, result.values, nz);
            }
            else
            {
                result.rowind = colind;
                result.colind = rowind;
                result.values = values;
            }

            return result;
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

        internal void Invalidate()
        {
            rowind = null;
            colind = null;
            values = null;
        }
    }
}
