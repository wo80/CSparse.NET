
namespace CSparse.Tests.Complex.Factorization
{
    using CSparse.Complex;
    using CSparse.Complex.Factorization;
    using Microsoft.VisualStudio.TestTools.UnitTesting;
    using System;
    using System.Numerics;

    [TestClass]
    public class SparseQRTest
    {
        private const double EPS = 1.0e-6;

        [TestMethod]
        public void TestSolve()
        {
            // Load matrix from a file.
            var A = ResourceLoader.Get<Complex>("general-40x40.mat");

            Assert.AreEqual(A.RowCount, A.ColumnCount);

            // Create test data.
            var x = Helper.CreateTestVector(A.ColumnCount);
            var b = Helper.Multiply(A, x);
            var r = Vector.Clone(b);

            var qr = new SparseQR(A, ColumnOrdering.MinimumDegreeAtA);

            qr.Solve(b, x);

            // Compute residual r = b - Ax.
            A.Multiply(-1.0, x, 1.0, r);

            Assert.IsTrue(Vector.Norm(r) < EPS);
        }

        [TestMethod]
        public void TestSolveOverdetermined()
        {
            // Load matrix from a file.
            var A = ResourceLoader.Get<Complex>("general-40x20.mat");

            int m = A.RowCount;
            int n = A.ColumnCount;

            Assert.IsTrue(m > n);

            // Create test data.
            var x = Helper.CreateTestVector(n);
            var b = Helper.Multiply(A, x);
            var r = Vector.Clone(b);

            var qr = new SparseQR(A, ColumnOrdering.MinimumDegreeAtA);

            // Compute min norm(Ax - b).
            qr.Solve(b, x);

            // Compute residual r = b - Ax.
            A.Multiply(-1.0, x, 1.0, r);

            Assert.IsTrue(Vector.Norm(r) < EPS);
        }

        [TestMethod]
        public void TestSolveUnderdetermined()
        {
            // Load matrix from a file.
            var A = ResourceLoader.Get<Complex>("general-20x40.mat");

            int m = A.RowCount;
            int n = A.ColumnCount;

            Assert.IsTrue(m < n);

            // Create test data.
            var x = Helper.CreateTestVector(n);
            var b = Helper.Multiply(A, x);
            var r = Vector.Clone(b);

            var qr = new SparseQR(A, ColumnOrdering.MinimumDegreeAtA);

            // Assuming A has full rank m, we have N(A) = n - m degrees of freedom.
            // Compute solution x with min norm(Ax - b).
            qr.Solve(b, x);

            // Compute residuals.
            A.Multiply(-1.0, x, 1.0, r);

            Assert.IsTrue(Vector.Norm(r) < EPS);
        }
    }
}
