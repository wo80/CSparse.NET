namespace CSparse.Tests.Complex.Factorization
{

    using System;
    using CSparse.Complex;
    using CSparse.Complex.Factorization;
    using NUnit.Framework;

    using Complex = System.Numerics.Complex;

    public class SparseLDLTest
    {
        private const double EPS = 1.0e-6;

        [Test]
        public void TestSolve()
        {
            // Load matrix from a file.
            var A = ResourceLoader.Get<Complex>("hermitian-40-spd.mat");

            // Create test data.
            var x = Helper.CreateTestVector(A.ColumnCount);
            var b = Helper.Multiply(A, x);
            var r = Vector.Clone(b);

            var ldl = SparseLDL.Create(A, ColumnOrdering.MinimumDegreeAtPlusA);

            // Solve Ax = b.
            ldl.Solve(b, x);

            // Compute residual r = b - Ax.
            A.Multiply(-1.0, x, 1.0, r);

            Assert.That(Vector.Norm(r.Length, r) < EPS, Is.True);
        }

        [Test]
        public void TestSolveNonSpd()
        {
            // Load matrix from a file.
            var A = ResourceLoader.Get<Complex>("hermitian-40.mat");

            // Create test data.
            var x = Helper.CreateTestVector(A.ColumnCount);
            var b = Helper.Multiply(A, x);
            var r = Vector.Clone(b);

            var ldl = SparseLDL.Create(A, ColumnOrdering.MinimumDegreeAtPlusA);

            // Solve Ax = b.
            ldl.Solve(b, x);

            // Compute residual r = b - Ax.
            A.Multiply(-1.0, x, 1.0, r);

            Assert.That(Vector.Norm(r.Length, r) < EPS, Is.True);
        }

        [Test]
        public void TestRefactorize()
        {
            var A = ResourceLoader.Get<Complex>("hermitian-40-spd.mat");

            var ldl = SparseLDL.Create(A, ColumnOrdering.MinimumDegreeAtPlusA);

            // Same sparsity pattern, different values (B = 2*A) : the numeric
            // refactorization must reuse the cached symbolic analysis.
            var B = A.Clone();
            var bv = B.Values;
            for (int i = 0; i < bv.Length; i++) bv[i] *= 2.0;

            ldl.Refactorize(B);

            var x = Helper.CreateTestVector(B.ColumnCount);
            var b = Helper.Multiply(B, x);
            var r = Vector.Clone(b);

            ldl.Solve(b, x);
            B.Multiply(-1.0, x, 1.0, r);

            Assert.That(Vector.Norm(r.Length, r) < EPS, Is.True);

            // Dimension guard.
            var small = new SparseMatrix(3, 3, 0);
            Assert.Throws<ArgumentException>(() => ldl.Refactorize(small));
        }
    }
}
