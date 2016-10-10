// -----------------------------------------------------------------------
// <copyright file="SparseMatrixStorage.cs">
// Copyright (c) 2012-2016, Christian Woltering
// </copyright>
// -----------------------------------------------------------------------

namespace CSparse.Storage
{
    using CSparse.Properties;
    using System;
    using System.Collections.Generic;

    /// <summary>
    /// Abstract base class for sparse matrix storage implementations.
    /// </summary>
    public abstract class SparseMatrixStorage<T> : ISparseMatrixStorage<T>
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

        protected readonly int nrows;
        protected readonly int ncols;

        /// <inheritdoc />
        public int RowCount
        {
            get { return nrows; }
        }

        /// <inheritdoc />
        public int ColumnCount
        {
            get { return ncols; }
        }

        /// <inheritdoc />
        public abstract int NonZerosCount
        {
            get;
        }

        /// <summary>
        /// Initializes a new instance of the SparseMatrixStorage class.
        /// </summary>
        protected SparseMatrixStorage(int rowCount, int columnCount)
        {
            if (rowCount < 0 || columnCount < 0)
            {
                throw new ArgumentOutOfRangeException(Resources.MatrixDimensionNonNegative);
            }

            this.nrows = rowCount;
            this.ncols = columnCount;
        }

        /// <inheritdoc />
        public abstract T At(int row, int column);

        /// <inheritdoc />
        public abstract void Clear();

        /// <inheritdoc />
        public abstract double Norm(int which);

        /// <inheritdoc />
        public abstract int Keep(Func<int, int, T, double, bool> func, double tolerance);

        /// <inheritdoc />
        public abstract IEnumerable<Tuple<int, int, T>> EnumerateIndexed();

        /// <summary>
        /// Removes numerically zero entries from a matrix.
        /// </summary>
        /// <param name="tolerance">Drop tolerance.</param>
        /// <returns>New number of non-zero entries in A.</returns>
        public abstract int DropZeros(double tolerance = 0.0);

        /// <summary>
        /// Computes the product y = alpha * A * x + beta * y.
        /// </summary>
        public abstract void Multiply(T alpha, T[] x, T beta, T[] y);

        /// <summary>
        /// Computes the product y = alpha * A^t * x + beta * y.
        /// </summary>
        public abstract void TransposeMultiply(T alpha, T[] x, T beta, T[] y);

        #region Storage equality

        /// <summary>
        /// Indicates whether the current object is equal to another object of the same type.
        /// </summary>
        /// <param name="other">
        /// An object to compare with this object.
        /// </param>
        /// <returns>
        /// <c>true</c> if the current object is equal to the <paramref name="other"/> parameter; otherwise, <c>false</c>.
        /// </returns>
        public virtual bool Equals(ISparseMatrixStorage<T> other)
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
        /// Check two matrices for equality.
        /// </summary>
        public abstract bool Equals(ISparseMatrixStorage<T> other, double tolerance);

        /// <summary>
        /// Determines whether the specified <see cref="T:System.Object"/> is equal to the current <see cref="T:System.Object"/>.
        /// </summary>
        /// <returns>
        /// true if the specified <see cref="T:System.Object"/> is equal to the current <see cref="T:System.Object"/>; otherwise, false.
        /// </returns>
        /// <param name="obj">The <see cref="T:System.Object"/> to compare with the current <see cref="T:System.Object"/>. </param>
        public override sealed bool Equals(object obj)
        {
            return Equals(obj as ISparseMatrixStorage<T>);
        }

        /// <summary>
        /// Serves as a hash function for a particular type.
        /// </summary>
        /// <returns>
        /// A hash code for the current <see cref="T:System.Object"/>.
        /// </returns>
        public abstract override int GetHashCode();

        #endregion
    }
}
