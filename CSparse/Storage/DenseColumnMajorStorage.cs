namespace CSparse.Storage
{
    using CSparse;
    using CSparse.Properties;
    using System;
    using System.Collections.Generic;
    using System.Numerics;

    /// <summary>
    /// Dense column-major matrix storage.
    /// </summary>
    /// <typeparam name="T">Supported data types are <c>double</c> and <see cref="T:Complex"/>.</typeparam>
    [Serializable]
    public abstract class DenseColumnMajorStorage<T> : Matrix<T>
        where T : struct, IEquatable<T>, IFormattable
    {
        /// <summary>
        /// Gets the numerical values in column-major order.
        /// </summary>
        public T[] Values;

        /// <summary>
        /// Return the matrix value at position (i, j).
        /// </summary>
        /// <param name="i">The row index.</param>
        /// <param name="j">The column index.</param>
        public T this[int i, int j]
        {
            get { return Values[(j * rowCount) + i]; }
            set { Values[(j * rowCount) + i] = value; }
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="DenseColumnMajorStorage{T}"/> class.
        /// </summary>
        /// <param name="rows">The number of rows.</param>
        /// <param name="columns">The number of columns.</param>
        /// <remarks>
        /// A new array for the matrix values will be allocated (size <c>rows * columns</c>).
        /// </remarks>
        public DenseColumnMajorStorage(int rows, int columns)
            : base(rows, columns)
        {
            this.Values = new T[rows * columns];
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="DenseColumnMajorStorage{T}"/> class.
        /// </summary>
        /// <param name="rows">The number of rows.</param>
        /// <param name="columns">The number of columns.</param>
        /// <param name="values">The values array (minimum size <c>rows * columns</c>)</param>
        public DenseColumnMajorStorage(int rows, int columns, T[] values)
            : base(rows, columns)
        {
            if (values.Length < rows * columns)
            {
                throw new ArgumentException("Invalid values array size.", nameof(values));
            }

            this.Values = values;
        }

        #region Public static functions

        /// <summary>
        /// Create a new dense matrix as a copy of the given other matrix.
        /// </summary>
        public static DenseColumnMajorStorage<T> OfMatrix(Matrix<T> matrix)
        {
            return OfIndexed(matrix.RowCount, matrix.ColumnCount, matrix.EnumerateIndexed());
        }

        /// <summary>
        /// Create a new dense matrix as a copy of the given two-dimensional array.
        /// </summary>
        public static DenseColumnMajorStorage<T> OfArray(T[,] array)
        {
            int rows = array.GetLength(0);
            int columns = array.GetLength(1);

            var A = Create(rows, columns);

            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    A.At(i, j, array[i, j]);
                }
            }

            return A;
        }

        /// <summary>
        /// Create a new dense matrix as a copy of the given indexed enumerable.
        /// </summary>
        public static DenseColumnMajorStorage<T> OfIndexed(int rows, int columns, IEnumerable<Tuple<int, int, T>> enumerable)
        {
            var A = Create(rows, columns);

            foreach (var item in enumerable)
            {
                A.At(item.Item1, item.Item2, item.Item3);
            }

            return A;
        }

        /// <summary>
        /// Create a new dense matrix as a copy of the given array (row-major order).
        /// </summary>
        public static DenseColumnMajorStorage<T> OfRowMajor(int rows, int columns, T[] rowMajor)
        {
            var A = Create(rows, columns);

            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    A.At(i, j, rowMajor[(i * columns) + j]);
                }
            }

            return A;
        }

        /// <summary>
        /// Create a new dense matrix as a copy of the given array (column-major order).
        /// </summary>
        public static DenseColumnMajorStorage<T> OfColumnMajor(int rows, int columns, T[] columnMajor)
        {
            var A = Create(rows, columns);

            Array.Copy(columnMajor, A.Values, rows * columns);

            return A;
        }

        /// <summary>
        /// Create a new square dense matrix with the diagonal as a copy of the given array.
        /// </summary>
        public static DenseColumnMajorStorage<T> OfDiagonalArray(T[] diagonal)
        {
            int order = diagonal.Length;

            var A = Create(order, order);
            
            for (int i = 0; i < order; i++)
            {
                A.At(i, i, diagonal[i]);
            }

            return A;
        }

        /// <summary>
        /// Create a new square dense matrix and initialize each diagonal value to the same provided value.
        /// </summary>
        public static DenseColumnMajorStorage<T> CreateDiagonal(int order, T value)
        {
            var A = Create(order, order);

            for (int i = 0; i < order; i++)
            {
                A.At(i, i, value);
            }

            return A;
        }

        /// <summary>
        /// Create a new square dense identity matrix where each diagonal value is set to one.
        /// </summary>
        public static DenseColumnMajorStorage<T> CreateIdentity(int order)
        {
            return CreateDiagonal(order, One);
        }

        #endregion

        /// <inheritdoc />
        public override T At(int row, int column)
        {
            return Values[(column * rowCount) + row];
        }

        /// <summary>
        /// Sets the element (without range checking).
        /// </summary>
        public void At(int row, int column, T value)
        {
            Values[(column * rowCount) + row] = value;
        }

        /// <inheritdoc />
        public override T[] Row(int row)
        {
            var target = new T[columnCount];

            Row(row, target);

            return target;
        }

        /// <inheritdoc />
        public override T[] Column(int column)
        {
            var target = new T[rowCount];

            Column(column, target);

            return target;
        }

        /// <inheritdoc />
        public override void Row(int row, T[] target)
        {
            for (int i = 0; i < columnCount; i++)
            {
                target[i] = Values[(i * rowCount) + row];
            }
        }

        /// <inheritdoc />
        public override void Column(int column, T[] target)
        {
            Array.Copy(Values, column * rowCount, target, 0, rowCount);
        }

        /// <summary>
        /// Copy values from array to matrix row.
        /// </summary>
        /// <param name="row">The row index.</param>
        /// <param name="values">The new values.</param>
        public void SetRow(int row, T[] values)
        {
            var target = this.Values;

            for (int i = 0; i < columnCount; i++)
            {
                target[(i * rowCount) + row] = values[i];
            }
        }

        /// <summary>
        /// Copy values from array to matrix column.
        /// </summary>
        /// <param name="column">The column index.</param>
        /// <param name="values">The new values.</param>
        public void SetColumn(int column, T[] values)
        {
            var target = this.Values;

            Array.Copy(values, 0, target, column * rowCount, rowCount);
        }

        #region Linear Algebra (Vector)

        #endregion

        #region Linear Algebra (Matrix)

        /// <summary>
        /// Returns the transpose of this matrix.
        /// </summary>
        public virtual DenseColumnMajorStorage<T> Transpose()
        {
            var result = DenseColumnMajorStorage<T>.Create(columnCount, rowCount);
            this.Transpose(result);
            return result;
        }

        /// <summary>
        /// Transpose this matrix and store the result in given matrix.
        /// </summary>
        public virtual void Transpose(DenseColumnMajorStorage<T> result)
        {
            var target = result.Values;

            for (int j = 0; j < ColumnCount; j++)
            {
                var index = j * RowCount;
                for (int i = 0; i < RowCount; i++)
                {
                    target[(i * ColumnCount) + j] = Values[index + i];
                }
            }
        }

        /// <summary>
        /// Adds two dense matrices, C = A + B.
        /// </summary>
        public DenseColumnMajorStorage<T> Add(DenseColumnMajorStorage<T> other)
        {
            int m = this.rowCount;
            int n = this.columnCount;

            // check inputs
            if (m != other.RowCount || n != other.ColumnCount)
            {
                throw new ArgumentException(Resources.MatrixDimensions, "other");
            }

            var result = DenseColumnMajorStorage<T>.Create(m, n);

            Add(other, result);

            return result;
        }

        /// <summary>
        /// Adds two dense matrices, C = A + B.
        /// </summary>
        /// <param name="other">The matrix added to this instance.</param>
        /// <param name="result">Contains the sum.</param>
        public abstract void Add(DenseColumnMajorStorage<T> other, DenseColumnMajorStorage<T> result);

        /// <summary>
        /// Dense matrix multiplication, C = A*B
        /// </summary>
        /// <param name="other">The dense matrix multiplied with this instance.</param>
        /// <returns>C = A*B</returns>
        public DenseColumnMajorStorage<T> Multiply(DenseColumnMajorStorage<T> other)
        {
            int m = this.rowCount;
            int n = other.columnCount;

            // check inputs
            if (this.columnCount != other.RowCount)
            {
                throw new ArgumentException(Resources.MatrixDimensions, "other");
            }

            var result = DenseColumnMajorStorage<T>.Create(m, n);

            Multiply(other, result);

            return result;
        }

        /// <summary>
        /// Dense matrix multiplication, C = A*B
        /// </summary>
        /// <param name="other">The dense matrix multiplied to this instance.</param>
        /// <param name="result">The product matrix.</param>
        public abstract void Multiply(DenseColumnMajorStorage<T> other, DenseColumnMajorStorage<T> result);

        /// <summary>
        /// Dense matrix multiplication, C = A*B
        /// </summary>
        /// <param name="other">The dense matrix multiplied to this instance.</param>
        /// <param name="options">Parallel options (optional).</param>
        /// <returns>C = A*B</returns>
        public DenseColumnMajorStorage<T> ParallelMultiply(DenseColumnMajorStorage<T> other, System.Threading.Tasks.ParallelOptions options = null)
        {
            int m = this.rowCount;
            int n = other.columnCount;

            // check inputs
            if (this.columnCount != other.RowCount)
            {
                throw new ArgumentException(Resources.MatrixDimensions, "other");
            }

            var result = DenseColumnMajorStorage<T>.Create(m, n);

            ParallelMultiply(other, result, options);

            return result;
        }

        /// <summary>
        /// Dense matrix multiplication, C = A*B
        /// </summary>
        /// <param name="other">The dense matrix multiplied to this instance.</param>
        /// <param name="result">The product matrix.</param>
        /// <param name="options">Parallel options (optional).</param>
        public virtual void ParallelMultiply(DenseColumnMajorStorage<T> other, DenseColumnMajorStorage<T> result, System.Threading.Tasks.ParallelOptions options = null)
        {
            Multiply(other, result);
        }

        /// <summary>
        /// Pointwise multiplies this matrix with another matrix and stores the result into the result matrix.
        /// </summary>
        /// <param name="other">The matrix to pointwise multiply with this one.</param>
        /// <param name="result">The matrix to store the result of the pointwise multiplication.</param>
        public abstract void PointwiseMultiply(DenseColumnMajorStorage<T> other, DenseColumnMajorStorage<T> result);

        #endregion

        /// <summary>
        /// Returns a clone of this matrix.
        /// </summary>
        public abstract DenseColumnMajorStorage<T> Clone();

        /// <summary>
        /// Returns a new matrix containing the upper triangle of this matrix.
        /// </summary>
        /// <returns>The upper triangle of this matrix.</returns>
        public virtual DenseColumnMajorStorage<T> UpperTriangle()
        {
            var result = Create(RowCount, ColumnCount);

            for (var row = 0; row < RowCount; row++)
            {
                for (var column = row; column < ColumnCount; column++)
                {
                    result.At(row, column, At(row, column));
                }
            }

            return result;
        }

        /// <summary>
        /// Returns a new matrix containing the lower triangle of this matrix.
        /// </summary>
        /// <returns>The lower triangle of this matrix.</returns>
        public virtual DenseColumnMajorStorage<T> LowerTriangle()
        {
            var result = Create(RowCount, ColumnCount);

            for (var row = 0; row < RowCount; row++)
            {
                for (var column = 0; column <= row && column < ColumnCount; column++)
                {
                    result.At(row, column, At(row, column));
                }
            }

            return result;
        }

        /// <summary>
        /// Puts the lower triangle of this matrix into the result matrix.
        /// </summary>
        /// <param name="result">Where to store the lower triangle.</param>
        /// <exception cref="ArgumentNullException">If <paramref name="result"/> is <see langword="null" />.</exception>
        /// <exception cref="ArgumentException">If the result matrix's dimensions are not the same as this matrix.</exception>
        public virtual void LowerTriangle(DenseColumnMajorStorage<T> result)
        {
            if (result == null)
            {
                throw new ArgumentNullException(nameof(result));
            }

            if (result.RowCount != RowCount || result.ColumnCount != ColumnCount)
            {
                throw new ArgumentException(Resources.MatrixDimensions, "result");
            }

            for (var row = 0; row < RowCount; row++)
            {
                for (var column = 0; column < ColumnCount; column++)
                {
                    result.At(row, column, row >= column ? At(row, column) : Zero);
                }
            }
        }

        /// <summary>
        /// Puts the upper triangle of this matrix into the result matrix.
        /// </summary>
        /// <param name="result">Where to store the lower triangle.</param>
        /// <exception cref="ArgumentNullException">If <paramref name="result"/> is <see langword="null" />.</exception>
        /// <exception cref="ArgumentException">If the result matrix's dimensions are not the same as this matrix.</exception>
        public virtual void UpperTriangle(DenseColumnMajorStorage<T> result)
        {
            if (result == null)
            {
                throw new ArgumentNullException(nameof(result));
            }

            if (result.RowCount != RowCount || result.ColumnCount != ColumnCount)
            {
                throw new ArgumentException(Resources.MatrixDimensions, "result");
            }

            for (var row = 0; row < RowCount; row++)
            {
                for (var column = 0; column < ColumnCount; column++)
                {
                    result.At(row, column, row <= column ? At(row, column) : Zero);
                }
            }
        }

        /// <summary>
        /// Returns a sub-matrix with values in given range.
        /// </summary>
        /// <param name="rowIndex">The row to start copying to.</param>
        /// <param name="rowCount">The number of rows to copy. Must be positive.</param>
        /// <param name="columnIndex">The column to start copying to.</param>
        /// <param name="columnCount">The number of columns to copy. Must be positive.</param>
        public DenseColumnMajorStorage<T> SubMatrix(int rowIndex, int rowCount, int columnIndex, int columnCount)
        {
            var result = DenseColumnMajorStorage<T>.Create(rowCount, columnCount);

            CopySubMatrixTo(result, rowIndex, 0, rowCount, columnIndex, 0, columnCount);

            return result;
        }

        /// <summary>
        /// Copies the values of a given matrix into a region in this matrix.
        /// </summary>
        /// <param name="rowIndex">The row to start copying to.</param>
        /// <param name="columnIndex">The column to start copying to.</param>
        /// <param name="subMatrix">The sub-matrix to copy from.</param>
        public void SetSubMatrix(int rowIndex, int columnIndex, DenseColumnMajorStorage<T> subMatrix)
        {
            subMatrix.CopySubMatrixTo(this, 0, rowIndex, subMatrix.RowCount, 0, columnIndex, subMatrix.ColumnCount);
        }

        /// <summary>
        /// Copies the values of a given matrix into a region in this matrix.
        /// </summary>
        /// <param name="rowIndex">The row to start copying to.</param>
        /// <param name="rowCount">The number of rows to copy. Must be positive.</param>
        /// <param name="columnIndex">The column to start copying to.</param>
        /// <param name="columnCount">The number of columns to copy. Must be positive.</param>
        /// <param name="subMatrix">The sub-matrix to copy from.</param>
        public void SetSubMatrix(int rowIndex, int rowCount, int columnIndex, int columnCount, DenseColumnMajorStorage<T> subMatrix)
        {
            subMatrix.CopySubMatrixTo(this, 0, rowIndex, rowCount, 0, columnIndex, columnCount);
        }

        /// <inheritdoc />
        public override void Clear()
        {
            Array.Clear(Values, 0, rowCount * columnCount);
        }

        /// <inheritdoc />
        public override IEnumerable<Tuple<int, int, T>> EnumerateIndexed()
        {
            for (int row = 0; row < rowCount; row++)
            {
                for (int column = 0; column < columnCount; column++)
                {
                    yield return new Tuple<int, int, T>(row, column, Values[(column * rowCount) + row]);
                }
            }
        }

        private void CopySubMatrixTo(DenseColumnMajorStorage<T> target,
            int sourceRowIndex, int targetRowIndex, int rowCount,
            int sourceColumnIndex, int targetColumnIndex, int columnCount)
        {
            if (target == null)
            {
                throw new ArgumentNullException(nameof(target));
            }

            if (rowCount == 0 || columnCount == 0)
            {
                return;
            }

            if (ReferenceEquals(this, target))
            {
                throw new NotSupportedException();
            }

            // TODO: Validate sub-matrix range.

            for (int j = sourceColumnIndex, jj = targetColumnIndex; j < sourceColumnIndex + columnCount; j++, jj++)
            {
                Array.Copy(Values, j * RowCount + sourceRowIndex, target.Values, jj * target.RowCount + targetRowIndex, rowCount);
            }
        }

        #region Internal methods

        internal static DenseColumnMajorStorage<T> Create(int rowCount, int columnCount)
        {
            if (typeof(T) == typeof(double))
            {
                return new CSparse.Double.DenseMatrix(rowCount, columnCount)
                    as DenseColumnMajorStorage<T>;
            }

            if (typeof(T) == typeof(Complex))
            {
                return new CSparse.Complex.DenseMatrix(rowCount, columnCount)
                    as DenseColumnMajorStorage<T>;
            }

            throw new NotSupportedException();
        }

        #endregion
    }
}
