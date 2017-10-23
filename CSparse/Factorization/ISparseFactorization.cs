// -----------------------------------------------------------------------
// <copyright file="ISparseFactorization.cs">
// Copyright (c) 2012-2016, Christian Woltering
// </copyright>
// -----------------------------------------------------------------------

namespace CSparse.Factorization
{
    using System;

    /// <summary>
    /// Interface for factorization methods.
    /// </summary>
    /// <typeparam name="T">Supported data types are <c>double</c> and <see cref="Complex"/>.</typeparam>
    public interface ISparseFactorization<T> : ISolver<T>
        where T : struct, IEquatable<T>, IFormattable
    {
        /// <summary>
        /// Gets the total number of non-zero entries (all factors).
        /// </summary>
        int NonZerosCount { get; }
    }
}
