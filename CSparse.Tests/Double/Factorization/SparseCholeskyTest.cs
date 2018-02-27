
namespace CSparse.Tests.Double.Factorization
{
    using CSparse.Double;
    using CSparse.Double.Factorization;
    using NUnit.Framework;
    using System;

    public class SparseCholeskyTest
    {
        private const double EPS = 1.0e-6;

        [Test]
        public void TestSolve()
        {
            // Load matrix from a file.
            var A = ResourceLoader.Get<double>("symmetric-40-spd.mat");

            // Create test data.
            var x = Helper.CreateTestVector(A.ColumnCount);
            var b = Helper.Multiply(A, x);
            var r = Vector.Clone(b);

            var chol = SparseCholesky.Create(A, ColumnOrdering.MinimumDegreeAtPlusA);

            // Solve Ax = b.
            chol.Solve(b, x);

            // Compute residual r = b - Ax.
            A.Multiply(-1.0, x, 1.0, r);

            Assert.IsTrue(Vector.Norm(r) < EPS);
        }

        [Test]
        public void TestConstructorThrowsOnNonSpd()
        {
            // Load matrix from a file.
            var A = ResourceLoader.Get<double>("symmetric-40.mat");

            Assert.Throws<Exception>(() =>
            {
                var chol = SparseCholesky.Create(A, ColumnOrdering.MinimumDegreeAtPlusA);
            });
        }

        [Test]
        public void TestEmptyFactorize()
        {
            var A = new SparseMatrix(0, 0, 0);

            var chol = SparseCholesky.Create(A, ColumnOrdering.MinimumDegreeAtPlusA);

            Assert.NotNull(chol);
            Assert.IsTrue(chol.NonZerosCount == 0);
        }
    }
}
