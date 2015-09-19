
namespace CSparse.Tests.Double
{
    using CSparse.Storage;

    static class Helper
    {
        public static double[] CreateTestVector(int n)
        {
            var x = new double[n];

            for (int i = 0; i < n; i++)
            {
                x[i] = 1.0 + i / (n - 1.0);
            }

            return x;
        }

        public static double[] Multiply(ISparseMatrixStorage<double> A, double[] x)
        {
            var b = new double[A.RowCount];

            A.Multiply(1.0, x, 0.0, b);

            return b;
        }
    }
}
