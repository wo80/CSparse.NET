
namespace CSparse.Tests.Complex
{
    using CSparse.Storage;
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

        public static Complex[] Multiply(ISparseMatrixStorage<Complex> A, Complex[] x)
        {
            var b = new Complex[A.RowCount];

            A.Multiply(1.0, x, 0.0, b);

            return b;
        }
    }
}
