namespace CSparse.Storage
{
    using CSparse.Properties;
    using System;
    using System.Collections.Generic;
    using System.Numerics;

    /// <summary>
    /// Compressed sparse column storage.
    /// </summary>
    /// <typeparam name="T"></typeparam>
    [Serializable]
    public abstract class CompressedColumnStorage<T> : Matrix<T>
        where T : struct, IEquatable<T>, IFormattable
    {
        /// <summary>
        /// Gets or sets a value indicating whether the storage should be
        /// automatically resized to non-zeros count. Defaults to true.
        /// </summary>
        public static bool AutoTrimStorage { get; set; } = true;

        /// <summary>
        /// Row pointers with last entry equal number of non-zeros (size = RowCount + 1)
        /// </summary>
        public int[] ColumnPointers;

        /// <summary>
        /// Column indices (size >= NonZerosCount)
        /// </summary>
        public int[] RowIndices;

        /// <summary>
        /// Numerical values (size >= NonZerosCount)
        /// </summary>
        public T[] Values;

        /// <summary>
        /// Gets the number of non-zero entries.
        /// </summary>
        public int NonZerosCount
        {
            get { return ColumnPointers[columns]; }
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="CompressedColumnStorage{T}"/> class.
        /// </summary>
        /// <param name="rowCount">The number of rows.</param>
        /// <param name="columnCount">The number of columns.</param>
        /// <remarks>By default, no arrays are allocated. The user will have to assign the storage arrays.</remarks>
        protected CompressedColumnStorage(int rowCount, int columnCount)
            : base(rowCount, columnCount)
        {
            // No array initialization here
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="CompressedColumnStorage{T}"/> class.
        /// </summary>
        /// <param name="rowCount">The number of rows.</param>
        /// <param name="columnCount">The number of columns.</param>
        /// <param name="valueCount">The number of non-zero values.</param>
        protected CompressedColumnStorage(int rowCount, int columnCount, int valueCount)
            : base(rowCount, columnCount)
        {
            if (valueCount < 0)
            {
                throw new ArgumentOutOfRangeException(nameof(valueCount), Resources.ValueNonNegative);
            }

            this.ColumnPointers = new int[columnCount + 1];
            this.RowIndices = new int[valueCount];
            this.Values = new T[valueCount];
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="CompressedColumnStorage{T}"/> class. Based on other CCS arrays
        /// </summary>
        /// <param name="rowCount">The number of rows.</param>
        /// <param name="columnCount">The number of columns.</param>
        /// <param name="values"></param>
        /// <param name="rowIndices"></param>
        /// <param name="columnPointers"></param>
        /// <remarks>The provided arrays will be re-used (not cloned).</remarks>
        public CompressedColumnStorage(int rowCount, int columnCount, T[] values, int[] rowIndices, int[] columnPointers)
            : base(rowCount, columnCount)
        {
            if (columnPointers.Length != columnCount + 1)
            {
                throw new ArgumentException("columnPointers.Length must equal (columnCount + 1)");
            }

            if (values.Length != rowIndices.Length)
            {
                throw new ArgumentException("rowIndices.Length must equal values.Length");
            }

            this.ColumnPointers = columnPointers;
            this.RowIndices = rowIndices;
            this.Values = values;
        }

        #region Public static functions

        /// <summary>
        /// Create a new sparse matrix as a copy of the given other matrix.
        /// </summary>
        public static CompressedColumnStorage<T> OfMatrix(Matrix<T> matrix)
        {
            var c = Converter.FromEnumerable_<T>(matrix.EnumerateIndexed(), matrix.RowCount, matrix.ColumnCount);

            return Converter.ToCompressedColumnStorage_(c);
        }

        /// <summary>
        /// Create a new sparse matrix as a copy of the given two-dimensional array.
        /// </summary>
        public static CompressedColumnStorage<T> OfArray(T[,] array)
        {
            var c = Converter.FromDenseArray_(array);

            return Converter.ToCompressedColumnStorage_(c);
        }

        /// <summary>
        /// Create a new sparse matrix as a copy of the given two-dimensional array.
        /// </summary>
        public static CompressedColumnStorage<T> OfJaggedArray(T[][] array)
        {
            var c = Converter.FromJaggedArray_(array);

            return Converter.ToCompressedColumnStorage_(c);
        }

        /// <summary>
        /// Create a new sparse matrix as a copy of the given coordinate storage.
        /// </summary>
        public static CompressedColumnStorage<T> OfIndexed(CoordinateStorage<T> coordinateStorage)
        {
            return Converter.ToCompressedColumnStorage_(coordinateStorage);
        }

        /// <summary>
        /// Create a new sparse matrix as a copy of the given indexed enumerable.
        /// </summary>
        public static CompressedColumnStorage<T> OfIndexed(int rows, int columns, IEnumerable<Tuple<int, int, T>> enumerable)
        {
            var c = Converter.FromEnumerable_<T>(enumerable, rows, columns);

            return Converter.ToCompressedColumnStorage_(c);
        }

        /// <summary>
        /// Create a new sparse matrix as a copy of the given array (row-major).
        /// </summary>
        public static CompressedColumnStorage<T> OfRowMajor(int rows, int columns, T[] rowMajor)
        {
            var c = Converter.FromRowMajorArray_<T>(rowMajor, rows, columns);

            return Converter.ToCompressedColumnStorage_(c);
        }

        /// <summary>
        /// Create a new sparse matrix as a copy of the given array (column-major).
        /// </summary>
        public static CompressedColumnStorage<T> OfColumnMajor(int rows, int columns, T[] columnMajor)
        {
            var c = Converter.FromColumnMajorArray_<T>(columnMajor, rows, columns);

            return Converter.ToCompressedColumnStorage_(c);
        }

        /// <summary>
        /// Create a new square sparse matrix with the diagonal as a copy of the given array.
        /// </summary>
        public static CompressedColumnStorage<T> OfDiagonalArray(T[] diagonal)
        {
            int order = diagonal.Length;

            var A = Create(order, order, order);

            var ap = A.ColumnPointers;
            var ai = A.RowIndices;
            var ax = A.Values;

            for (int i = 0; i < order; i++)
            {
                ap[i] = i;
                ai[i] = i;
                ax[i] = diagonal[i];
            }

            ap[order] = order;

            return A;
        }

        /// <summary>
        /// Create a new square sparse matrix and initialize each diagonal value to the same provided value.
        /// </summary>
        public static CompressedColumnStorage<T> CreateDiagonal(int order, T value)
        {
            var A = Create(order, order, order);

            var ap = A.ColumnPointers;
            var ai = A.RowIndices;
            var ax = A.Values;

            for (int i = 0; i < order; i++)
            {
                ap[i] = i;
                ai[i] = i;
                ax[i] = value;
            }

            ap[order] = order;

            return A;
        }

        /// <summary>
        /// Create a new square sparse identity matrix where each diagonal value is set to One.
        /// </summary>
        public static CompressedColumnStorage<T> CreateIdentity(int order)
        {
            return CreateDiagonal(order, One);
        }

        #endregion

        /// <inheritdoc />
        public override T At(int row, int column)
        {
            int index = ColumnPointers[column];
            int length = ColumnPointers[column + 1] - index;
            int pos = Array.BinarySearch(RowIndices, index, length, row);
            return pos >= 0 ? Values[pos] : Zero;
        }

        /// <inheritdoc />
        public override void Clear()
        {
            Array.Clear(ColumnPointers, 0, ColumnPointers.Length);
            Array.Clear(Values, 0, Values.Length);
        }

        /// <inheritdoc />
        public override T[] Row(int rowIndex)
        {
            var target = new T[ColumnCount];

            Row(rowIndex, target);

            return target;
        }

        /// <inheritdoc />
        public override void Row(int rowIndex, T[] target)
        {
            if (target.Length != ColumnCount)
            {
                throw new Exception();
            }

            int i, n = ColumnCount;

            var ap = ColumnPointers;
            var ai = RowIndices;
            var ax = Values;

            for (int k = 0; k < n; k++)
            {
                // Check if columns contain row index.
                i = Array.BinarySearch(ai, ap[k], ap[k + 1] - ap[k], rowIndex);

                if (i >= 0)
                {
                    target[k] = ax[i];
                }
            }
        }

        /// <inheritdoc />
        public override T[] Column(int columnIndex)
        {
            var target = new T[RowCount];

            Column(columnIndex, target);

            return target;
        }

        /// <inheritdoc />
        public override void Column(int columnIndex, T[] target)
        {
            if (target.Length != RowCount)
            {
                throw new Exception();
            }

            var ap = ColumnPointers;
            var ai = RowIndices;
            var ax = Values;

            int colEnd = ap[columnIndex + 1];

            for (int k = ap[columnIndex]; k < colEnd; k++)
            {
                target[ai[k]] = ax[k];
            }
        }

        #region Linear Algebra (Matrix)

        /// <summary>
        /// Returns the transpose of this matrix.
        /// </summary>
        public CompressedColumnStorage<T> Transpose()
        {
            return this.Transpose(false);
        }

        /// <summary>
        /// Transpose this matrix and store the result in given matrix.
        /// </summary>
        /// <param name="result">Storage for the tranposed matrix.</param>
        public void Transpose(CompressedColumnStorage<T> result)
        {
            this.Transpose(result, false);
        }

        /// <summary>
        /// Returns the transpose of this matrix.
        /// </summary>
        /// <param name="storage">A value indicating, whether the transpose should be done on storage level (without complex conjugation).</param>
        public CompressedColumnStorage<T> Transpose(bool storage)
        {
            var result = CompressedColumnStorage<T>.Create(columns, rows, this.NonZerosCount);
            this.Transpose(result, storage);
            return result;
        }

        /// <summary>
        /// Transpose this matrix and store the result in given matrix.
        /// </summary>
        /// <param name="result">Storage for the tranposed matrix.</param>
        /// <param name="storage">A value indicating, whether the transpose should be done on storage level (without complex conjugation).</param>
        public virtual void Transpose(CompressedColumnStorage<T> result, bool storage)
        {
            int i, j, p;

            var cx = result.Values;
            var cp = result.ColumnPointers;
            var ci = result.RowIndices;

            int[] w = new int[rows];

            for (p = 0; p < ColumnPointers[columns]; p++)
            {
                // Row counts.
                w[RowIndices[p]]++;
            }

            // Row pointers.
            Helper.CumulativeSum(cp, w, rows);

            for (i = 0; i < columns; i++)
            {
                for (p = ColumnPointers[i]; p < ColumnPointers[i + 1]; p++)
                {
                    j = w[RowIndices[p]]++;

                    // Place A(i,j) as entry C(j,i)
                    ci[j] = i;
                    cx[j] = Values[p];
                }
            }
        }

        /// <summary>
        /// Adds two matrices in CSC format, C = A + B, where A is current instance.
        /// </summary>
        public CompressedColumnStorage<T> Add(CompressedColumnStorage<T> other)
        {
            int m = this.rows;
            int n = this.columns;

            // check inputs
            if (m != other.RowCount || n != other.ColumnCount)
            {
                throw new ArgumentException(Resources.MatrixDimensions);
            }

            var result = CompressedColumnStorage<T>.Create(m, n, this.NonZerosCount + other.NonZerosCount);

            var one = Helper.OneOf<T>();

            Add(one, one, other, result);

            return result;
        }

        /// <summary>
        /// Adds two matrices, C = alpha*A + beta*B, where A is current instance.
        /// </summary>
        /// <param name="alpha">Scalar factor for A, current instance.</param>
        /// <param name="beta">Scalar factor for B, other instance.</param>
        /// <param name="other">The matrix added to this instance.</param>
        /// <param name="result">Contains the sum.</param>
        /// <remarks>
        /// The <paramref name="result"/> matrix has to be fully initialized and provide enough
        /// space for the nonzero entries of the sum. An upper bound is the sum of the nonzeros
        /// count of A and B.
        /// </remarks>
        public abstract void Add(T alpha, T beta, CompressedColumnStorage<T> other,
            CompressedColumnStorage<T> result);

        /// <summary>
        /// Sparse matrix multiplication, C = A*B
        /// </summary>
        /// <param name="other">The sparse matrix multiplied to this instance.</param>
        /// <returns>C = A*B</returns>
        public CompressedColumnStorage<T> Multiply(CompressedColumnStorage<T> other)
        {
            int m = this.rows;
            int n = other.ColumnCount;

            var result = CompressedColumnStorage<T>.Create(m, n, this.NonZerosCount + other.NonZerosCount);

            Multiply(other, result);

            return result;
        }

        /// <summary>
        /// Sparse matrix multiplication, C = A*B
        /// </summary>
        /// <param name="other">The sparse matrix multiplied to this instance.</param>
        /// <param name="result">Contains the matrix product.</param>
        /// <remarks>
        /// The <paramref name="result"/> matrix has to be fully initialized, but doesn't have
        /// to provide enough space for the nonzero entries of the product. The storage will be
        /// automatically expanded if necessary.
        /// </remarks>
        public abstract void Multiply(CompressedColumnStorage<T> other, CompressedColumnStorage<T> result);

        /// <summary>
        /// Sparse matrix multiplication, C = A*B
        /// </summary>
        /// <param name="other">The sparse matrix multiplied to this instance.</param>
        /// <param name="options">Parallel options (optional).</param>
        /// <returns>C = A*B</returns>
        public virtual CompressedColumnStorage<T> ParallelMultiply(CompressedColumnStorage<T> other, System.Threading.Tasks.ParallelOptions options = null)
        {
            return Multiply(other);
        }

        #endregion

        /// <summary>
        /// Filter matrix values.
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
        public abstract int Keep(Func<int, int, T, bool> func);

        /// <summary>
        /// Removes numerically zero entries from a matrix.
        /// </summary>
        /// <param name="tolerance">Drop tolerance (default is 0.0)</param>
        /// <returns>The new number of nonzero entries.</returns>
        public abstract int DropZeros(double tolerance = 0.0);

        /// <summary>
        /// Returns a clone of this matrix.
        /// </summary>
        /// <param name="values">If true (default), the values are copied.</param>
        public CompressedColumnStorage<T> Clone(bool values = true)
        {
            int rows = this.RowCount;
            int columns = this.ColumnCount;

            int nnz = this.NonZerosCount;

            var ap = this.ColumnPointers;
            var ai = this.RowIndices;

            var result = Create(rows, columns, values ? nnz : 0);

            if (values)
            {
                Array.Copy(this.Values, 0, result.Values, 0, nnz);
            }
            else if (nnz > 0)
            {
                // Fix size of row indices array in case values == false.
                result.RowIndices = new int[nnz];
            }

            Buffer.BlockCopy(ap, 0, result.ColumnPointers, 0, (columns + 1) * Constants.SizeOfInt);
            Buffer.BlockCopy(ai, 0, result.RowIndices, 0, nnz * Constants.SizeOfInt);

            return result;
        }

        /// <inheritdoc />
        public override IEnumerable<Tuple<int, int, T>> EnumerateIndexed()
        {
            int n = this.ColumnCount;

            var ax = this.Values;
            var ap = this.ColumnPointers;
            var ai = this.RowIndices;

            for (int i = 0; i < n; i++)
            {
                var end = ap[i + 1];
                for (var j = ap[i]; j < end; j++)
                {
                    yield return new Tuple<int, int, T>(ai[j], i, ax[j]);
                }
            }
        }

        /// <summary>
        /// Evaluates whether this matrix is symmetric.
        /// </summary>
        public virtual bool IsSymmetric()
        {
            int n = this.ColumnCount;

            if (this.RowCount != n)
            {
                return false;
            }

            var ax = this.Values;
            var ap = this.ColumnPointers;
            var ai = this.RowIndices;

            // If we assume that columns are sorted, the symmetry check can be
            // made more efficient, checking only entries above the diagonal.

            for (var i = 0; i < n; i++)
            {
                int end = ap[i + 1];

                for (var j = ap[i]; j < end; j++)
                {
                    if (!ax[j].Equals(At(i, ai[j])))
                    {
                        return false;
                    }
                }
            }

            return true;
        }

        /// <summary>
        /// Permute the rows of the matrix.
        /// </summary>
        /// <param name="perm">Permutation matrix P.</param>
        /// <param name="target">The target storage (must be fully initialized to match the source storage).</param>
        public void PermuteRows(int[] perm, CompressedColumnStorage<T> target)
        {
            var bx = target.Values;
            var bp = target.ColumnPointers;
            var bi = target.RowIndices;

            if (target.rows != rows || target.columns != columns)
            {
                throw new ArgumentException(Resources.InvalidDimensions, "target");
            }

            if (perm.Length != rows)
            {
                throw new ArgumentException("Invalid permutation length.", "perm");
            }

            PermuteRows(Values, ColumnPointers, RowIndices, bx, bp, bi, perm);

            Helper.SortIndices(target);
        }

        /// <summary>
        /// Permute the rows of the matrix.
        /// </summary>
        /// <param name="perm">Permutation matrix P.</param>
        public void PermuteRows(int[] perm)
        {
            PermuteRows(perm, this);
        }

        /// <summary>
        /// Permute the columns of the matrix.
        /// </summary>
        /// <param name="perm">Permutation matrix P.</param>
        /// <param name="target">The target storage (must be fully initialized to match the source storage).</param>
        public void PermuteColumns(int[] perm, CompressedColumnStorage<T> target)
        {
            var bx = target.Values;
            var bp = target.ColumnPointers;
            var bi = target.RowIndices;

            if (ReferenceEquals(this, target))
            {
                throw new ArgumentException("Cannot use this instance as target.", "target");
            }

            if (target.rows != rows || target.columns != columns)
            {
                throw new ArgumentException(Resources.InvalidDimensions, "target");
            }

            if (perm.Length != columns)
            {
                throw new ArgumentException("Invalid permutation length.", "perm");
            }

            PermuteColumns(Values, ColumnPointers, RowIndices, bx, bp, bi, perm);

            Helper.SortIndices(target);
        }

        /// <summary>
        /// Permute the columns of the matrix.
        /// </summary>
        /// <param name="perm">Permutation matrix P.</param>
        public CompressedColumnStorage<T> PermuteColumns(int[] perm)
        {
            var result = Create(RowCount, columns, Values.Length);

            PermuteColumns(perm, result);

            return result;
        }


        /// <summary>
        /// Returns the positions of the diagonal elements of a sparse matrix.
        /// </summary>
        /// <param name="throwOnMissingDiag"></param>
        /// <returns></returns>
        public int[] FindDiagonalIndices(bool throwOnMissingDiag = false)
        {
            int n = this.ColumnCount;

            var ap = this.ColumnPointers;
            var ai = this.RowIndices;

            int[] diag = new int[n];

            for (int i = 0; i < n; i++)
            {
                diag[i] = Array.BinarySearch(ai, ap[i], ap[i + 1] - ap[i], i);

                if (diag[i] < 0 && throwOnMissingDiag)
                {
                    throw new Exception("Missing diagonal entry on row " + (i + 1));
                }
            }

            return diag;
        }

        #region Permutation methods

        /// <summary>
        /// Permutes the columns of a matrix in CSC format, B = A * P, where P represents
        /// a permutation matrix.
        /// </summary>
        /// <param name="ax">Input matrix values.</param>
        /// <param name="ai">Input matrix row pointers.</param>
        /// <param name="aj">Input matrix column indices.</param>
        /// <param name="bx">Output matrix values.</param>
        /// <param name="bi">Output matrix row pointers.</param>
        /// <param name="bj">Output matrix column indices.</param>
        /// <param name="perm">Permutation array of length ColumnCount.</param>
        /// <remarks>
        /// The permutation P is defined through the array perm: for each j,
        /// perm(j) represents the destination row number of row number j:
        /// 
        /// a(i,j) in the original matrix becomes a(perm(i),j) in the output matrix.
        /// </remarks>
        protected void PermuteColumns(T[] ax, int[] ai, int[] aj, T[] bx, int[] bi, int[] bj, int[] perm)
        {
            int k;

            // Determine pointers for output matix. 
            for (int i = 0; i < columns; i++)
            {
                k = perm[i];
                bi[k + 1] = ai[i + 1] - ai[i];
            }

            // Get pointers from lengths
            bi[0] = 0;
            for (int i = 0; i < columns; i++)
            {
                bi[i + 1] += bi[i];
            }

            // Copying
            for (int i = 0; i < columns; i++)
            {
                // Old row = i, new row = perm(i), k = new pointer
                k = bi[perm[i]];
                for (int j = ai[i]; j < ai[i + 1]; j++)
                {
                    bj[k] = aj[j];
                    bx[k] = ax[j];
                    k = k + 1;
                }
            }
        }

        /// <summary>
        /// Permute the rows of a matrix in CSC format, B = P * A, where P represents
        /// a permutation matrix. 
        /// </summary>
        /// <param name="ax">Input matrix values.</param>
        /// <param name="ai">Input matrix row pointers.</param>
        /// <param name="aj">Input matrix column indices.</param>
        /// <param name="bx">Output matrix values.</param>
        /// <param name="bi">Output matrix row pointers.</param>
        /// <param name="bj">Output matrix column indices.</param>
        /// <param name="perm">Permutation array of length RowCount.</param>
        /// <param name="copy">Copy matrix values (not needed if used 'in place').</param>
        /// <remarks>
        /// The permutation matrix P maps column j into column perm(j), i.e., 
        /// on return a(i,j) in the original matrix becomes a(i,perm(j)) in the
        /// output matrix.
        /// 
        /// Notes:
        /// 
        /// 1. This routine is in place: aj, bj can be the same.
        /// 2. If the matrix is initially sorted (by increasing column number) 
        ///    then bx, bi, bj may not be on return.
        /// </remarks>
        protected void PermuteRows(T[] ax, int[] ai, int[] aj, T[] bx, int[] bi, int[] bj,
            int[] perm, bool copy = false)
        {
            int i, nnz = ai[columns];

            for (i = 0; i < nnz; i++)
            {
                bj[i] = perm[aj[i]];
            }

            if (copy)
            {
                Array.Copy(ax, bx, nnz);
                Array.Copy(ai, bi, columns);
            }
        }

        #endregion

        #region Internal methods

        internal static CompressedColumnStorage<T> Create(int rowCount, int columnCount)
        {
            if (typeof(T) == typeof(double))
            {
                return new CSparse.Double.SparseMatrix(rowCount, columnCount)
                    as CompressedColumnStorage<T>;
            }

            if (typeof(T) == typeof(Complex))
            {
                return new CSparse.Complex.SparseMatrix(rowCount, columnCount)
                    as CompressedColumnStorage<T>;
            }

            throw new NotSupportedException();
        }

        internal static CompressedColumnStorage<T> Create(int rowCount, int columnCount, int valueCount)
        {
            if (typeof(T) == typeof(double))
            {
                return new CSparse.Double.SparseMatrix(rowCount, columnCount, valueCount)
                    as CompressedColumnStorage<T>;
            }

            if (typeof(T) == typeof(Complex))
            {
                return new CSparse.Complex.SparseMatrix(rowCount, columnCount, valueCount)
                    as CompressedColumnStorage<T>;
            }

            throw new NotSupportedException();
        }

        /// <summary>
        /// Change the max # of entries sparse matrix
        /// </summary>
        /// <param name="size"></param>
        /// <returns></returns>
        internal bool Resize(int size)
        {
            if (size <= 0)
            {
                size = this.ColumnPointers[columns];
            }

            // TODO: check available memory
            //       and throw OutOfMemoryException before trying to resize.

            Array.Resize<int>(ref this.RowIndices, size);
            Array.Resize<T>(ref this.Values, size);

            return true;
        }

        internal abstract void Cleanup();

        internal abstract int Scatter(int j, T beta, int[] w, T[] x, int mark, CompressedColumnStorage<T> mat, int nzz);

        #endregion

        #region Storage equality

        /// <summary>
        /// Serves as a hash function for a particular type.
        /// </summary>
        /// <returns>
        /// A hash code for the current <see cref="T:System.Object"/>.
        /// </returns>
        public override int GetHashCode()
        {
            var hashNum = Math.Min(this.NonZerosCount, 25);
            int hash = 17;
            int i, p, k = 0;
            unchecked
            {
                for (i = 0; i < columns; i++)
                {
                    for (p = ColumnPointers[i]; p < ColumnPointers[i + 1]; p++)
                    {
                        hash = hash * 31 + Values[p].GetHashCode();

                        if (++k > hashNum)
                        {
                            return hash;
                        }
                    }
                }
            }
            return hash;
        }

        #endregion
    }
}
