// -----------------------------------------------------------------------
// <copyright file="Helper.cs">
// Copyright (c) 2006-2016, Timothy A. Davis
// Copyright (c) 2012-2016, Christian Woltering
// </copyright>
// -----------------------------------------------------------------------

namespace CSparse
{
    using System;

    internal static class Helper
    {
        /// <summary>
        /// Cumulative sum of given array.
        /// </summary>
        /// <param name="sum">Output: cumulative sum of counts</param>
        /// <param name="counts">input array, overwritten with sum</param>
        /// <param name="size">length of counts</param>
        /// <returns>sum[size] (non-zeros)</returns>
        public static int CumulativeSum(int[] sum, int[] counts, int size)
        {
            int i, nz = 0;

            for (i = 0; i < size; i++)
            {
                sum[i] = nz;
                nz += counts[i];
                counts[i] = sum[i]; // also copy p[0..n-1] back into c[0..n-1]
            }

            sum[size] = nz;

            return nz;
        }

        #region Generic code helper

        /// <summary>
        /// Sets the value of <c>1.0</c> for type T.
        /// </summary>
        /// <typeparam name="T">The type to return the value of 1.0 of.</typeparam>
        /// <returns>The value of <c>1.0</c> for type T.</returns>
        public static T OneOf<T>()
        {
            if (typeof(T) == typeof(double))
            {
                return (T)(object)1.0d;
            }

            if (typeof(T) == typeof(System.Numerics.Complex))
            {
                return (T)(object)System.Numerics.Complex.One;
            }

            throw new NotSupportedException();
        }

        /// <summary>
        /// Sets the value of <c>0.0</c> for type T.
        /// </summary>
        /// <typeparam name="T">The type to return the value of 0.0 of.</typeparam>
        /// <returns>The value of <c>0.0</c> for type T.</returns>
        public static T ZeroOf<T>()
        {
            if (typeof(T) == typeof(double))
            {
                return (T)(object)0.0d;
            }

            if (typeof(T) == typeof(System.Numerics.Complex))
            {
                return (T)(object)System.Numerics.Complex.Zero;
            }

            throw new NotSupportedException();
        }

        #endregion
    }
}
