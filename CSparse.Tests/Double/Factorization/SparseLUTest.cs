namespace CSparse.Tests.Double.Factorization
{
    using CSparse.Double;
    using CSparse.Double.Factorization;
    using CSparse.Storage;
    using NUnit.Framework;
    using System;

    public class SparseLUTest
    {
        private const double EPS = 1.0e-6;

        [Test]
        public void TestSolve()
        {
            // Load matrix from a file.
            var A = ResourceLoader.Get<double>("general-40x40.mat");

            // Create test data.
            var x = Helper.CreateTestVector(A.ColumnCount);
            var b = Helper.Multiply(A, x);
            var r = Vector.Clone(b);

            // Create LU factorization.
            var lu = SparseLU.Create(A, ColumnOrdering.MinimumDegreeAtPlusA, 1.0);

            // Solve Ax = b.
            lu.Solve(b, x);

            // Compute residual r = b - Ax.
            A.Multiply(-1.0, x, 1.0, r);

            Assert.That(Vector.Norm(r.Length, r) < EPS, Is.True);

            // Test exceptions:

            var e1 = Assert.Throws<ArgumentNullException>(() => lu.Solve(b, null));
            var e2 = Assert.Throws<ArgumentNullException>(() => lu.Solve(null, x));

            Assert.That(e1.ParamName, Is.EqualTo("result"));
            Assert.That(e2.ParamName, Is.EqualTo("input"));
        }

        [Test]
        public void TestSolveTranspose()
        {
            // Load matrix from a file.
            var A = ResourceLoader.Get<double>("general-40x40.mat");

            var AT = A.Transpose();

            // Create test data.
            var x = Helper.CreateTestVector(A.ColumnCount);
            var b = Helper.Multiply(AT, x);
            var r = Vector.Clone(b);

            // Create LU factorization.
            var lu = SparseLU.Create(A, ColumnOrdering.MinimumDegreeAtPlusA, 1.0);

            // Solve A'x = b.
            lu.SolveTranspose(b, x);

            // Compute residual r = b - A'x.
            AT.Multiply(-1.0, x, 1.0, r);

            Assert.That(Vector.Norm(r.Length, r) < EPS, Is.True);
        }

        [Test]
        public void TestEmptyFactorize()
        {
            var A = new SparseMatrix(0, 0, 0);

            var lu = SparseLU.Create(A, ColumnOrdering.MinimumDegreeAtPlusA, 1.0);

            Assert.That(lu, Is.Not.Null);
            Assert.That(lu.NonZerosCount == 0, Is.True);
        }

        [Test]
        public void TestRefactorize()
        {
            // Load matrix from a file.
            var A = ResourceLoader.Get<double>("general-40x40.mat");

            // Symbolic + numeric factorization of A.
            var lu = SparseLU.Create(A, ColumnOrdering.MinimumDegreeAtPlusA, 1.0);

            // Same sparsity pattern, different values (B = 2*A) : numeric refactorization
            // must reuse the cached symbolic ordering and still solve correctly.
            var B = A.Clone();
            var bv = B.Values;
            for (int i = 0; i < bv.Length; i++) bv[i] *= 2.0;

            lu.Refactorize(B, 1.0);

            var x = Helper.CreateTestVector(B.ColumnCount);
            var b = Helper.Multiply(B, x);
            var r = Vector.Clone(b);

            // Solve Bx = b with the refactorized decomposition.
            lu.Solve(b, x);

            // Residual r = b - Bx.
            B.Multiply(-1.0, x, 1.0, r);

            Assert.That(Vector.Norm(r.Length, r) < EPS, Is.True);

            // Dimension guard.
            var small = new SparseMatrix(3, 3, 0);
            Assert.Throws<ArgumentException>(() => lu.Refactorize(small, 1.0));
        }

        [Test]
        public void TestRefactorizeNoTrim()
        {
            // With AutoTrimStorage disabled, the L/U buffers are kept across
            // re-factorizations (no Resize(0) trim) — the result must stay correct.
            var trim = CompressedColumnStorage<double>.AutoTrimStorage;
            try
            {
                CompressedColumnStorage<double>.AutoTrimStorage = false;

                var A = ResourceLoader.Get<double>("general-40x40.mat");
                var lu = SparseLU.Create(A, ColumnOrdering.MinimumDegreeAtPlusA, 1.0);

                var B = A.Clone();
                var bv = B.Values;
                for (int i = 0; i < bv.Length; i++) bv[i] *= 2.0;

                lu.Refactorize(B, 1.0);

                var x = Helper.CreateTestVector(B.ColumnCount);
                var b = Helper.Multiply(B, x);
                var r = Vector.Clone(b);

                lu.Solve(b, x);
                B.Multiply(-1.0, x, 1.0, r);

                Assert.That(Vector.Norm(r.Length, r) < EPS, Is.True);
            }
            finally
            {
                CompressedColumnStorage<double>.AutoTrimStorage = trim;
            }
        }
    }
}
