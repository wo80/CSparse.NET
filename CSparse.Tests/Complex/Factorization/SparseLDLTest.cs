namespace CSparse.Tests.Complex.Factorization
{

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

            Assert.IsTrue(Vector.Norm(r.Length, r) < EPS);
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

            Assert.IsTrue(Vector.Norm(r.Length, r) < EPS);
        }
    }
}
