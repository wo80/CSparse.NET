namespace CSparse.Double
{
    using System;

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
        public static void Copy(double[] src, double[] dst)
        {
            Buffer.BlockCopy(src, 0, dst, 0, src.Length * Constants.SizeOfDouble);
        }

        /// <summary>
        /// Copy one vector to another.
        /// </summary>
        /// <param name="n">Number of values to copy.</param>
        /// <param name="src">The source array.</param>
        /// <param name="dst">The destination array.</param>
        public static void Copy(int n, double[] src, double[] dst)
        {
            Buffer.BlockCopy(src, 0, dst, 0, n * Constants.SizeOfDouble);
        }

        /// <summary>
        /// Create a new vector.
        /// </summary>
        public static double[] Create(int length, double value)
        {
            double[] result = new double[length];

            for (int i = 0; i < length; i++)
            {
                result[i] = value;
            }

            return result;
        }

        /// <summary>
        /// Clone the given vector.
        /// </summary>
        public static double[] Clone(double[] src)
        {
            double[] result = new double[src.Length];

            Buffer.BlockCopy(src, 0, result, 0, src.Length * Constants.SizeOfDouble);

            return result;
        }

        /// <summary>
        /// Set vector values to zero.
        /// </summary>
        public static void Clear(double[] x)
        {
            Array.Clear(x, 0, x.Length);
        }

        /// <summary>
        /// Computes the dot product of two vectors.
        /// </summary>
        [Obsolete("Use DotProduct(n, x, y).")]
        public static double DotProduct(double[] x, double[] y)
        {
            return DotProduct(x.Length, x, y);
        }

        /// <summary>
        /// Computes the dot product of two vectors.
        /// </summary>
        public static double DotProduct(int n, double[] x, double[] y)
        {
            double result = 0.0;

            for (int i = 0; i < n; i++)
            {
                result += x[i] * y[i];
            }

            return result;
        }

        /// <summary>
        /// Computes the pointwise product of two vectors.
        /// </summary>
        [Obsolete("Use PointwiseMultiply(n, x, y, target).")]
        public static void PointwiseMultiply(double[] x, double[] y, double[] target)
        {
            PointwiseMultiply(x.Length, x, y, target);
        }

        /// <summary>
        /// Computes the pointwise product of two vectors.
        /// </summary>
        public static void PointwiseMultiply(int n, double[] x, double[] y, double[] target)
        {
            for (int i = 0; i < n; i++)
            {
                target[i] = x[i] * y[i];
            }
        }

        /// <summary>
        /// Computes the norm of a vector, sqrt( x' * x ).
        /// </summary>
        [Obsolete("Use Norm(n, x).")]
        public static double Norm(double[] x)
        {
            return Norm(x.Length, x);
        }

        /// <summary>
        /// Computes the norm of a vector, sqrt( x' * x ).
        /// </summary>
        public static double Norm(int n, double[] x)
        {
            double result = 0.0;

            for (int i = 0; i < n; ++i)
            {
                result += x[i] * x[i];
            }

            return Math.Sqrt(result);
        }

        /// <summary>
        /// Computes the norm of a vector avoiding overflow, sqrt( x' * x ).
        /// </summary>
        [Obsolete("Use NormRobust(n, x).")]
        public static double NormRobust(double[] x)
        {
            return NormRobust(x.Length, x);
        }

        /// <summary>
        /// Computes the norm of a vector avoiding overflow, sqrt( x' * x ).
        /// </summary>
        public static double NormRobust(int n, double[] x)
        {
            double scale = 0.0, ssq = 1.0;

            for (int i = 0; i < n; ++i)
            {
                if (x[i] != 0.0)
                {
                    double absxi = Math.Abs(x[i]);
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
        public static void Scale(double a, double[] x)
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
        public static void Scale(int n, double a, double[] x, double[] target)
        {
            for (int i = 0; i < n; i++)
            {
                target[i] = a * x[i];
            }
        }

        /// <summary>
        /// Add a scaled vector to another vector, y = a * x + y.
        /// </summary>
        public static void Axpy(double a, double[] x, double[] y)
        {
            int length = x.Length;

            for (int i = 0; i < length; i++)
            {
                y[i] += a * x[i];
            }
        }

        /// <summary>
        /// Add two vectors, z = a * x + b * y.
        /// </summary>
        public static void Add(double a, double[] x, double b, double[] y, double[] target)
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
        public static void Add(int n, double a, double[] x, double[] y, double[] target)
        {
            for (int i = 0; i < n; i++)
            {
                target[i] = a * x[i] + y[i];
            }
        }

        /// <summary>
        /// Add two vectors, target = a * x + b * y.
        /// </summary>
        public static void Add(int n, double a, double[] x, double b, double[] y, double[] target)
        {
            for (int i = 0; i < n; i++)
            {
                target[i] = a * x[i] + b * y[i];
            }
        }

        /// <summary>
        /// Add three vectors, target = a * x + b * y + z.
        /// </summary>
        public static void Add(int n, double a, double[] x, double b, double[] y, double[] z, double[] target)
        {
            for (int i = 0; i < n; i++)
            {
                target[i] = a * x[i] + b * y[i] + z[i];
            }
        }

        /// <summary>
        /// Add three vectors, target = a * x + b * y + c * z.
        /// </summary>
        public static void Add(int n, double a, double[] x, double b, double[] y, double c, double[] z, double[] target)
        {
            for (int i = 0; i < n; i++)
            {
                target[i] = a * x[i] + b * y[i] + c * z[i];
            }
        }
    }
}
