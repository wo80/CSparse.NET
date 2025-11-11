namespace CSparse.Tests.Complex.Factorization
{
    using CSparse.Complex;
    using CSparse.Complex.Factorization;
    using NUnit.Framework;
    using System;
    using Complex = System.Numerics.Complex;

    public class SparseLUTest
    {
        private const double EPS = 1.0e-6;

        [Test]
        public void TestSolve()
        {
            // Load matrix from a file.
            var A = ResourceLoader.Get<Complex>("general-40x40.mat");

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
            var A = ResourceLoader.Get<Complex>("general-40x40.mat");

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
    }
}
