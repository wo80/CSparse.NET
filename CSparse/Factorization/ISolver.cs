namespace CSparse.Factorization
{
    using System;

    /// <summary>
    /// Classes that solve a system of linear equations, <c>Ax = b</c>.
    /// </summary>
    /// <typeparam name="T">Supported data types are <c>double</c> and <see cref="System.Numerics.Complex"/>.</typeparam>
    public interface ISolver<T> where T : struct, IEquatable<T>, IFormattable
    {
        /// <summary>
        /// Solves a system of linear equations, Ax = b.
        /// </summary>
        /// <param name="input">Right hand side b</param>
        /// <param name="result">Solution vector x.</param>
        void Solve(T[] input, T[] result);

        /// <summary>
        /// Solves a system of linear equations, Ax = b.
        /// </summary>
        /// <param name="input">Right hand side b</param>
        /// <param name="result">Solution vector x.</param>
        void Solve(ReadOnlySpan<T> input, Span<T> result);
    }
}
