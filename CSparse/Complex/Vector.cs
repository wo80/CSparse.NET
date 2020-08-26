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
        /// <param name="n">Number of values to copy.</param>
        /// <param name="src">The source array.</param>
        /// <param name="dst">The destination array.</param>
        public static void Copy(int n, Complex[] src, Complex[] dst)
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
        [Obsolete("Use DotProduct(n, x, y).")]
        public static Complex DotProduct(Complex[] x, Complex[] y)
        {
            return DotProduct(x.Length, x, y);
        }

        /// <summary>
        /// Computes the dot product of two vectors.
        /// </summary>
        public static Complex DotProduct(int n, Complex[] x, Complex[] y)
        {
            Complex result = 0.0;

            for (int i = 0; i < n; i++)
            {
                result += Complex.Conjugate(x[i]) * y[i];
            }

            return result;
        }

        /// <summary>
        /// Computes the pointwise product of two vectors.
        /// </summary>
        [Obsolete("Use PointwiseMultiply(n, x, y, target).")]
        public static void PointwiseMultiply(Complex[] x, Complex[] y, Complex[] target)
        {
            PointwiseMultiply(x.Length, x, y, target);
        }

        /// <summary>
        /// Computes the pointwise product of two vectors.
        /// </summary>
        public static void PointwiseMultiply(int n, Complex[] x, Complex[] y, Complex[] z)
        {
            for (int i = 0; i < n; i++)
            {
                z[i] = x[i] * y[i];
            }
        }

        /// <summary>
        /// Computes the norm of a vector, sqrt( x' * x ).
        /// </summary>
        [Obsolete("Use Norm(n, x).")]
        public static double Norm(Complex[] x)
        {
            return Norm(x.Length, x);
        }

        /// <summary>
        /// Computes the norm of a vector.
        /// </summary>
        public static double Norm(int n, Complex[] x)
        {
            double re, im, sum = 0.0;

            for (int i = 0; i < n; i++)
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
        [Obsolete("Use NormRobust(n, x).")]
        public static double NormRobust(Complex[] x)
        {
            return NormRobust(x.Length, x);
        }

        /// <summary>
        /// Computes the norm of a vector avoiding overflow, sqrt( x' * x ).
        /// </summary>
        public static double NormRobust(int n, Complex[] x)
        {
            double scale = 0.0, ssq = 1.0;

            for (int i = 0; i < n; ++i)
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
        /// Scales a vector by a given factor, x = a * x.
        /// </summary>
        public static void Scale(Complex a, Complex[] x)
        {
            int length = x.Length;

            for (int i = 0; i < length; i++)
            {
                x[i] *= a;
            }
        }

        /// <summary>
        /// Scales a vector by a given factor, target = a * x.
        /// </summary>
        public static void Scale(int n, Complex a, Complex[] x, Complex[] target)
        {
            for (int i = 0; i < n; i++)
            {
                target[i] = a * x[i];
            }
        }

        /// <summary>
        /// Add a scaled vector to another vector, y = a * x + y.
        /// </summary>
        public static void Axpy(Complex a, Complex[] x, Complex[] y)
        {
            int length = x.Length;

            for (int i = 0; i < length; i++)
            {
                y[i] += a * x[i];
            }
        }

        /// <summary>
        /// Add two scaled vectors, z = a * x + b * y.
        /// </summary>
        public static void Add(Complex a, Complex[] x, Complex b, Complex[] y, Complex[] target)
        {
            int length = x.Length;

            for (int i = 0; i < length; i++)
            {
                target[i] = a * x[i] + b * y[i];
            }
        }

        /// <summary>
        /// Add two vectors, target = a * x + y.
        /// </summary>
        public static void Add(int n, Complex a, Complex[] x, Complex[] y, Complex[] target)
        {
            for (int i = 0; i < n; i++)
            {
                target[i] = a * x[i] + y[i];
            }
        }

        /// <summary>
        /// Add two vectors, target = a * x + b * y.
        /// </summary>
        public static void Add(int n, Complex a, Complex[] x, Complex b, Complex[] y, Complex[] target)
        {
            for (int i = 0; i < n; i++)
            {
                target[i] = a * x[i] + b * y[i];
            }
        }

        /// <summary>
        /// Add three vectors, target = a * x + b * y + z.
        /// </summary>
        public static void Add(int n, Complex a, Complex[] x, Complex b, Complex[] y, Complex[] z, Complex[] target)
        {
            for (int i = 0; i < n; i++)
            {
                target[i] = a * x[i] + b * y[i] + z[i];
            }
        }

        /// <summary>
        /// Add three vectors, target = a * x + b * y + c * z.
        /// </summary>
        public static void Add(int n, Complex a, Complex[] x, Complex b, Complex[] y, Complex c, Complex[] z, Complex[] target)
        {
            for (int i = 0; i < n; i++)
            {
                target[i] = a * x[i] + b * y[i] + c * z[i];
            }
        }
    }
}
