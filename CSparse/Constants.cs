// -----------------------------------------------------------------------
// <copyright file="Constants.cs">
// Copyright (c) 2012-2015, Christian Woltering
// </copyright>
// -----------------------------------------------------------------------

namespace CSparse
{
    using System.Runtime.InteropServices;

    /// <summary>
    /// Constants used in the library.
    /// </summary>
    public static class Constants
    {
        /// <summary>
        /// The size of an int in bytes.
        /// </summary>
        public const int SizeOfInt = sizeof(int);

        /// <summary>
        /// The size of a double in bytes.
        /// </summary>
        public const int SizeOfDouble = sizeof(double);

        /// <summary>
        /// The default threshold used for matrix values comparision.
        /// </summary>
        public const double EqualsThreshold = 1e-12;

        /// <summary>
        /// The size of a Complex in bytes (should be 2 * SizeOfDouble).
        /// </summary>
        public static readonly int SizeOfComplex = Marshal.SizeOf(typeof(System.Numerics.Complex));
    }
}
