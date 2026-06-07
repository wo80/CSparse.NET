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

            Assert.That(Vector.Norm(r.Length, r) < EPS, Is.True);
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

            Assert.That(chol, Is.Not.Null);
            Assert.That(chol.NonZerosCount == 0, Is.True);
        }

        [Test]
        public void TestRefactorize()
        {
            var A = ResourceLoader.Get<double>("symmetric-40-spd.mat");

            var chol = SparseCholesky.Create(A, ColumnOrdering.MinimumDegreeAtPlusA);

            // Same sparsity pattern, different values (B = 2*A stays SPD) : the numeric
            // refactorization must reuse the cached symbolic analysis and still solve.
            var B = A.Clone();
            var bv = B.Values;
            for (int i = 0; i < bv.Length; i++) bv[i] *= 2.0;

            chol.Refactorize(B);

            var x = Helper.CreateTestVector(B.ColumnCount);
            var b = Helper.Multiply(B, x);
            var r = Vector.Clone(b);

            chol.Solve(b, x);
            B.Multiply(-1.0, x, 1.0, r);

            Assert.That(Vector.Norm(r.Length, r) < EPS, Is.True);

            // Dimension guard.
            var small = new SparseMatrix(3, 3, 0);
            Assert.Throws<ArgumentException>(() => chol.Refactorize(small));
        }
    }
}
