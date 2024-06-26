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
        /// The default threshold used for matrix values comparison.
        /// </summary>
        public const double EqualsThreshold = 1e-12;


        /// <summary>
        /// Machine epsilon for double precision.
        /// </summary>
        public const double MachineEpsilon = 2.2204460492503131e-16;

        /// <summary>
        /// The size of a Complex in bytes (should be 2 * SizeOfDouble).
        /// </summary>
        public static readonly int SizeOfComplex = Marshal.SizeOf(typeof(System.Numerics.Complex));
    }
}
