
namespace CSparse.Tests.Complex.Factorization
{
    using CSparse.Complex;
    using CSparse.Complex.Factorization;
    using NUnit.Framework;
    using System.Numerics;

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

            Assert.IsTrue(Vector.Norm(r) < EPS);
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

            Assert.IsTrue(Vector.Norm(r) < EPS);
        }

        [Test]
        public void TestEmptyFactorize()
        {
            var A = MatrixHelper.Load(0, 0);

            var lu = SparseLU.Create(A, ColumnOrdering.MinimumDegreeAtPlusA, 1.0);

            Assert.NotNull(lu);
            Assert.IsTrue(lu.NonZerosCount == 0);
        }
    }
}
