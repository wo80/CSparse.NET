namespace CSparse
{
    using System;

    /// <summary>
    /// Linear operator interface.
    /// </summary>
    /// <typeparam name="T">Supported data types are <c>double</c> and <see cref="System.Numerics.Complex"/>.</typeparam>
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
        /// Multiplies a (m-by-n) matrix by a vector, y = A*x. 
        /// </summary>
        /// <param name="x">Vector of length n (column count).</param>
        /// <param name="y">Vector of length m (row count), containing the result.</param>
        void Multiply(T[] x, T[] y);

        /// <summary>
        /// Multiplies a (m-by-n) matrix by a vector, y = A*x. 
        /// </summary>
        /// <param name="x">Vector of length n (column count).</param>
        /// <param name="y">Vector of length m (row count), containing the result.</param>
        void Multiply(ReadOnlySpan<T> x, Span<T> y);

        /// <summary>
        /// Multiplies a (m-by-n) matrix by a vector, y = alpha * A * x + beta * y.
        /// </summary>
        /// <param name="alpha">Scaling factor for vector x.</param>
        /// <param name="x">Vector of length n (column count).</param>
        /// <param name="beta">Scaling factor for vector y.</param>
        /// <param name="y">Vector of length m (row count), containing the result.</param>
        /// <remarks>
        /// Input values of vector <paramref name="y"/> will be accumulated.
        /// </remarks>
        void Multiply(T alpha, T[] x, T beta, T[] y);

        /// <summary>
        /// Multiplies a (m-by-n) matrix by a vector, y = alpha * A * x + beta * y.
        /// </summary>
        /// <param name="alpha">Scaling factor for vector x.</param>
        /// <param name="x">Vector of length n (column count).</param>
        /// <param name="beta">Scaling factor for vector y.</param>
        /// <param name="y">Vector of length m (row count), containing the result.</param>
        /// <remarks>
        /// Input values of vector <paramref name="y"/> will be accumulated.
        /// </remarks>
        void Multiply(T alpha, ReadOnlySpan<T> x, T beta, Span<T> y);

        /// <summary>
        /// Multiplies the transpose of a (m-by-n) matrix by a vector, y = A'*x. 
        /// </summary>
        /// <param name="x">Vector of length m (column count of A').</param>
        /// <param name="y">Vector of length n (row count of A'), containing the result.</param>
        void TransposeMultiply(T[] x, T[] y);

        /// <summary>
        /// Multiplies the transpose of a (m-by-n) matrix by a vector, y = A'*x. 
        /// </summary>
        /// <param name="x">Vector of length m (column count of A').</param>
        /// <param name="y">Vector of length n (row count of A'), containing the result.</param>
        void TransposeMultiply(ReadOnlySpan<T> x, Span<T> y);

        /// <summary>
        /// Multiplies the transpose of a (m-by-n) matrix by a vector, y = alpha * A^t * x + beta * y.
        /// </summary>
        /// <param name="alpha">Scaling factor for vector x.</param>
        /// <param name="x">Vector of length m (column count of A').</param>
        /// <param name="beta">Scaling factor for vector y.</param>
        /// <param name="y">Vector of length n (row count of A'), containing the result.</param>
        /// <remarks>
        /// Input values of vector <paramref name="y"/> will be accumulated.
        /// </remarks>
        void TransposeMultiply(T alpha, T[] x, T beta, T[] y);

        /// <summary>
        /// Multiplies the transpose of a (m-by-n) matrix by a vector, y = alpha * A^t * x + beta * y.
        /// </summary>
        /// <param name="alpha">Scaling factor for vector x.</param>
        /// <param name="x">Vector of length m (column count of A').</param>
        /// <param name="beta">Scaling factor for vector y.</param>
        /// <param name="y">Vector of length n (row count of A'), containing the result.</param>
        /// <remarks>
        /// Input values of vector <paramref name="y"/> will be accumulated.
        /// </remarks>
        void TransposeMultiply(T alpha, ReadOnlySpan<T> x, T beta, Span<T> y);
    }
}
