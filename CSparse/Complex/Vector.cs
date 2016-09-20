// -----------------------------------------------------------------------
// <copyright file="Vector.cs">
// Copyright (c) 2012-2016, Christian Woltering
// </copyright>
// -----------------------------------------------------------------------

namespace CSparse.Complex
{
    using System;
    using System.Numerics;

    /// <summary>
    /// Vector helper methods.
    /// </summary>
    public static class Vector
    {
        /// <summary>
        /// Copy one vector to another.
        /// </summary>
        /// <param name="src">The source array.</param>
        /// <param name="dst">The destination array.</param>
        public static void Copy(Complex[] src, Complex[] dst)
        {
            Array.Copy(src, dst, src.Length);
        }

        /// <summary>
        /// Copy one vector to another.
        /// </summary>
        /// <param name="src">The source array.</param>
        /// <param name="dst">The destination array.</param>
        /// <param name="n">Number of values to copy.</param>
        public static void Copy(Complex[] src, Complex[] dst, int n)
        {
            Array.Copy(src, dst, n);
        }

        /// <summary>
        /// Create a new vector.
        /// </summary>
        public static Complex[] Create(int length, Complex value)
        {
            Complex[] result = new Complex[length];

            for (int i = 0; i < length; i++)
            {
                result[i] = value;
            }

            return result;
        }

        /// <summary>
        /// Clone the given vector.
        /// </summary>
        public static Complex[] Clone(Complex[] src)
        {
            Complex[] result = new Complex[src.Length];

            Array.Copy(src, result, src.Length);

            return result;
        }
        
        /// <summary>
        /// Set vector values to zero.
        /// </summary>
        public static void Clear(Complex[] x)
        {
            Array.Clear(x, 0, x.Length);
        }

        /// <summary>
        /// Computes the dot product of two vectors.
        /// </summary>
        public static Complex DotProduct(Complex[] x, Complex[] y)
        {
            int length = x.Length;

            Complex result = 0.0;

            for (int i = 0; i < length; i++)
            {
                result += Complex.Conjugate(x[i]) * y[i];
            }

            return result;
        }

        /// <summary>
        /// Computes the norm of a vector.
        /// </summary>
        public static double Norm(Complex[] x)
        {
            double re, im, sum = 0.0;

            int length = x.Length;

            for (int i = 0; i < length; i++)
            {
                re = x[i].Real;
                im = x[i].Imaginary;
                sum += (re * re + im * im);
            }

            return Math.Sqrt(sum);
        }

        /// <summary>
        /// Computes the norm of a vector avoiding overflow, sqrt( x' * x ).
        /// </summary>
        public static double NormRobust(Complex[] x)
        {
            int length = x.Length;

            double scale = 0.0, ssq = 1.0;

            for (int i = 0; i < length; ++i)
            {
                if (x[i] != 0.0)
                {
                    double absxi = Complex.Abs(x[i]);
                    if (scale < absxi)
                    {
                        ssq = 1.0 + ssq * (scale / absxi) * (scale / absxi);
                        scale = absxi;
                    }
                    else
                    {
                        ssq += (absxi / scale) * (absxi / scale);
                    }
                }
            }

            return scale * Math.Sqrt(ssq);
        }

        /// <summary>
        /// Scales a vector by a given factor, x = alpha * x.
        /// </summary>
        public static void Scale(Complex alpha, Complex[] x)
        {
            int length = x.Length;

            for (int i = 0; i < length; i++)
            {
                x[i] *= alpha;
            }
        }

        /// <summary>
        /// Add a scaled vector t o another vector, y = y + alpha * x.
        /// </summary>
        public static void Axpy(Complex alpha, Complex[] x, Complex[] y)
        {
            int length = x.Length;

            for (int i = 0; i < length; i++)
            {
                y[i] += alpha * x[i];
            }
        }

        /// <summary>
        /// Add two scaled vectors, z = alpha * x + beta * y.
        /// </summary>
        public static void Add(Complex alpha, Complex[] x, Complex beta, Complex[] y, Complex[] z)
        {
            int length = x.Length;

            for (int i = 0; i < length; i++)
            {
                z[i] = alpha * x[i] + beta * y[i];
            }
        }
    }
}
