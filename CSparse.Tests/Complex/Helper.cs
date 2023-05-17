namespace CSparse.Tests.Complex
{
    using CSparse.Storage;
    using NUnit.Framework;
    using System.Numerics;

    static class Helper
    {
        public static Complex[] CreateTestVector(int n)
        {
            var x = new Complex[n];

            for (int i = 0; i < n; i++)
            {
                x[i] = 1.0 + i / (n - 1.0);
            }

            return x;
        }

        public static Complex[] Multiply(Matrix<Complex> A, Complex[] x)
        {
            var b = new Complex[A.RowCount];

            A.Multiply(1.0, x, 0.0, b);

            return b;
        }

        public static Complex[] TransposeMultiply(Matrix<Complex> A, Complex[] x, bool conjugate = true)
        {
            if (conjugate)
            {
                var conjA = ((CompressedColumnStorage<Complex>)A).Clone();
                var values = conjA.Values;
                for (int i = 0; i < values.Length; i++)
                {
                    values[i] = Complex.Conjugate(values[i]);
                }
                A = conjA;
            }

            var b = new Complex[A.RowCount];

            A.TransposeMultiply(1.0, x, 0.0, b);

            return b;
        }

        public static void AssertEqual(int length, Complex[] x, Complex[] y)
        {
            for (int i = 0; i < length; i++)
            {
                Assert.AreEqual(x[i].Real, y[i].Real);
                Assert.AreEqual(x[i].Imaginary, y[i].Imaginary);
            }
        }
    }
}
