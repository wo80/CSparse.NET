// -----------------------------------------------------------------------
// <copyright file="ILinearOperator.cs">
// Copyright (c) 2012-2016, Christian Woltering
// </copyright>
// -----------------------------------------------------------------------

namespace CSparse
{
    using System;

    /// <summary>
    /// Linear operator interface.
    /// </summary>
    /// <typeparam name="T"></typeparam>
    public interface ILinearOperator<T> where T : struct, IEquatable<T>, IFormattable
    {
        /// <summary>
        /// Gets the number of rows.
        /// </summary>
        int RowCount { get; }

        /// <summary>
        /// Gets the number of columns.
        /// </summary>
        int ColumnCount { get; }

        /// <summary>
        /// Computes the matrix-vector product y = alpha * A * x + beta * y.
        /// </summary>
        void Multiply(T alpha, T[] x, T beta, T[] y);

        /// <summary>
        /// Computes the matrix-vector product y = alpha * A^t * x + beta * y.
        /// </summary>
        void TransposeMultiply(T alpha, T[] x, T beta, T[] y);
    }
}
