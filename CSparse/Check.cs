
namespace CSparse
{
    using CSparse.Properties;
    using CSparse.Storage;
    using System;

    static class Check
    {
        public static void NotNull(object obj, string name)
        {
            if (obj == null)
            {
                throw new ArgumentNullException(name);
            }
        }

        public static void NotNaN(double value, string name)
        {
            if (double.IsNaN(value))
            {
                throw new ArgumentException(Resources.ValueNotNaN, name);
            }
        }

        public static void SquareMatrix<T>(CompressedColumnStorage<T> A, string name)
            where T : struct, IEquatable<T>, IFormattable
        {
            if (A.RowCount != A.ColumnCount)
            {
                throw new ArgumentException(Resources.MatrixSquare, name);
            }
        }

        public static void Permutation(int[] p, int n, string name)
        {
            if (p.Length < n)
            {
                throw new ArgumentException(Resources.InvalidDimensions, name);
            }

            if (!CSparse.Permutation.IsValid(p))
            {
                throw new ArgumentException(Resources.InvalidPermutation, name);
            }
        }

        public static void Dimension(int actual, int expected, string name)
        {
            if (actual != expected)
            {
                throw new ArgumentException(Resources.InvalidDimensions, name);
            }
        }
    }
}
