using CSparse.Properties;
using System;
using System.Collections.Generic;

[assembly: CLSCompliant(true)]

namespace CSparse
{
    /// <summary>
    /// Abstract base class for matrix implementations.
    /// </summary>
    /// <typeparam name="T">Supported data types are <c>double</c> and <see cref="T:Complex"/>.</typeparam>
    [Serializable]
    public abstract class Matrix<T> : ILinearOperator<T>
        where T : struct, IEquatable<T>, IFormattable
    {
        /// <summary>
        /// Zero value for T.
        /// </summary>
        protected static readonly T Zero = Helper.ZeroOf<T>();

        /// <summary>
        /// One value for T.
        /// </summary>
        protected static readonly T One = Helper.OneOf<T>();

        protected readonly int rows;
        protected readonly int columns;

        /// <inheritdoc />
        public int RowCount
        {
            get { return rows; }
        }

        /// <inheritdoc />
        public int ColumnCount
        {
            get { return columns; }
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="Matrix{T}"/> class.
        /// </summary>
        protected Matrix(int rowCount, int columnCount)
        {
            if (rowCount < 0)
            {
                throw new ArgumentOutOfRangeException(nameof(rowCount), Resources.MatrixDimensionNonNegative);
            }

            if (columnCount < 0)
            {
                throw new ArgumentOutOfRangeException(nameof(columnCount), Resources.MatrixDimensionNonNegative);
            }

            this.rows = rowCount;
            this.columns = columnCount;
        }

        /// <summary>
        /// Return the matrix value at position (row, column).
        /// </summary>
        /// <param name="row">The row index.</param>
        /// <param name="column">The column index.</param>
        /// <returns>Matrix value</returns>
        public abstract T At(int row, int column);

        /// <summary>
        /// Clears all values form the matrix.
        /// </summary>
        /// <remarks>
        /// The method does not release memory.
        /// </remarks>
        public abstract void Clear();

        /// <summary>
        /// Extract row from matrix.
        /// </summary>
        /// <param name="rowIndex">The column index to extract.</param>
        public abstract T[] Row(int rowIndex);

        /// <summary>
        /// Extract row from matrix.
        /// </summary>
        /// <param name="rowIndex">The column index to extract.</param>
        /// <param name="target">Dense array.</param>
        public abstract void Row(int rowIndex, T[] target);

        /// <summary>
        /// Extract column from matrix.
        /// </summary>
        /// <param name="columnIndex">The column index to extract.</param>
        public abstract T[] Column(int columnIndex);

        /// <summary>
        /// Extract column from matrix.
        /// </summary>
        /// <param name="columnIndex">The column index to extract.</param>
        /// <param name="target">Dense array.</param>
        public abstract void Column(int columnIndex, T[] target);

        /// <summary>
        /// Calculates the induced L1 norm of this matrix.
        /// </summary>
        /// <returns>The maximum absolute column sum of the matrix.</returns>
        public abstract double L1Norm();

        /// <summary>
        /// Calculates the induced infinity norm of this matrix.
        /// </summary>
        /// <returns>The maximum absolute row sum of the matrix.</returns>
        public abstract double InfinityNorm();

        /// <summary>
        /// Calculates the entry-wise Frobenius norm of this matrix.
        /// </summary>
        /// <returns>The square root of the sum of the squared values.</returns>
        public abstract double FrobeniusNorm();

        /// <summary>
        /// Enumerates all values of the matrix.
        /// </summary>
        /// <returns>Enumeration of tuples (i, j, a[i, j]).</returns>
        public abstract IEnumerable<Tuple<int, int, T>> EnumerateIndexed();

        /// <inheritdoc />
        public abstract void Multiply(T[] x, T[] y);

        /// <inheritdoc />
        public abstract void Multiply(T alpha, T[] x, T beta, T[] y);

        /// <inheritdoc />
        public abstract void TransposeMultiply(T[] x, T[] y);

        /// <inheritdoc />
        public abstract void TransposeMultiply(T alpha, T[] x, T beta, T[] y);

        #region Storage equality

        /// <summary>
        /// Indicates whether the current <see cref="Matrix{T}"/> object is equal to another matrix of the same type.
        /// </summary>
        /// <param name="other">A <see cref="Matrix{T}"/> object to compare with this matrix.</param>
        /// <returns>
        /// Returns <c>true</c> if the current object is equal to the <paramref name="other"/> parameter; otherwise, <c>false</c>.
        /// </returns>
        public virtual bool Equals(Matrix<T> other)
        {
            // Reject equality when the argument is null or has a different shape.
            if (other == null)
            {
                return false;
            }

            // Accept if the argument is the same object as this.
            if (ReferenceEquals(this, other))
            {
                return true;
            }

            return Equals(other, Constants.EqualsThreshold);
        }

        /// <summary>
        /// Indicates whether the current <see cref="Matrix{T}"/> object is equal to another matrix of the same type.
        /// </summary>
        /// <param name="other">A <see cref="Matrix{T}"/> object to compare with this matrix.</param>
        /// <param name="tolerance">Tolerance value for two matrix entries to be the same.</param>
        public abstract bool Equals(Matrix<T> other, double tolerance);

        /// <summary>
        /// Determines whether the specified <see cref="T:System.Object"/> is equal to the current <see cref="T:System.Object"/>.
        /// </summary>
        /// <returns>
        /// true if the specified <see cref="T:System.Object"/> is equal to the current <see cref="T:System.Object"/>; otherwise, false.
        /// </returns>
        /// <param name="obj">The <see cref="T:System.Object"/> to compare with the current <see cref="T:System.Object"/>. </param>
        public override sealed bool Equals(object obj)
        {
            return Equals(obj as Matrix<T>);
        }

        /// <summary>
        /// Serves as a hash function for a particular type.
        /// </summary>
        /// <returns>
        /// A hash code for the current <see cref="T:System.Object"/>.
        /// </returns>
        public override int GetHashCode()
        {
            var hashNum = Math.Min(RowCount * ColumnCount, 25);
            int hash = 17;
            unchecked
            {
                for (var i = 0; i < hashNum; i++)
                {
                    var col = i % ColumnCount;
                    var row = (i - col) / RowCount;
                    hash = hash * 31 + At(row, col).GetHashCode();
                }
            }
            return hash;
        }

        #endregion
    }
}
