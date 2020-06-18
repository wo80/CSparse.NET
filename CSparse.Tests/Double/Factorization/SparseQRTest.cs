namespace CSparse.Tests.Double.Factorization
{
    using CSparse.Double;
    using CSparse.Double.Factorization;
    using NUnit.Framework;

    public class SparseQRTest
    {
        private const double EPS = 1.0e-6;

        [Test]
        public void TestSolve()
        {
            // Load matrix from a file.
            var A = ResourceLoader.Get<double>("general-40x40.mat");

            Assert.AreEqual(A.RowCount, A.ColumnCount);

            // Create test data.
            var x = Helper.CreateTestVector(A.ColumnCount);
            var b = Helper.Multiply(A, x);
            var r = Vector.Clone(b);

            var qr = SparseQR.Create(A, ColumnOrdering.MinimumDegreeAtA);

            // Solve Ax = b.
            qr.Solve(b, x);

            // Compute residual r = b - Ax.
            A.Multiply(-1.0, x, 1.0, r);

            Assert.IsTrue(Vector.Norm(r) < EPS);
        }

        [Test]
        public void TestSolveTranspose()
        {
            // Load matrix from a file.
            var A = ResourceLoader.Get<double>("general-40x40.mat");

            Assert.AreEqual(A.RowCount, A.ColumnCount);

            var AT = A.Transpose();

            // Create test data.
            var x = Helper.CreateTestVector(A.ColumnCount);
            var b = Helper.Multiply(AT, x);
            var r = Vector.Clone(b);

            // Create LU factorization.
            var qr = SparseQR.Create(A, ColumnOrdering.MinimumDegreeAtA);

            // Solve A'x = b.
            qr.SolveTranspose(b, x);

            // Compute residual r = b - A'x.
            AT.Multiply(-1.0, x, 1.0, r);

            Assert.IsTrue(Vector.Norm(r) < EPS);
        }

        [Test]
        public void TestSolveOverdetermined()
        {
            // Load matrix from a file.
            var A = ResourceLoader.Get<double>("general-40x20.mat");

            Assert.IsTrue(A.RowCount > A.ColumnCount);

            // Create test data.
            var x = Helper.CreateTestVector(A.ColumnCount);
            var b = Helper.Multiply(A, x);
            var r = Vector.Clone(b);

            var qr = SparseQR.Create(A, ColumnOrdering.MinimumDegreeAtA);

            // Compute min norm(Ax - b).
            qr.Solve(b, x);

            // Compute residual r = b - Ax.
            A.Multiply(-1.0, x, 1.0, r);

            Assert.IsTrue(Vector.Norm(r) < EPS);
        }

        [Test]
        public void TestSolveTransposeOverdetermined()
        {
            // Load matrix from a file.
            var A = ResourceLoader.Get<double>("general-40x20.mat");

            Assert.IsTrue(A.RowCount > A.ColumnCount);

            var AT = A.Transpose();

            // Create test data.
            var x = Helper.CreateTestVector(A.RowCount);
            var b = Helper.Multiply(AT, x);
            var r = Vector.Clone(b);

            var qr = SparseQR.Create(A, ColumnOrdering.MinimumDegreeAtA);

            // Compute min norm(A'x - b).
            qr.SolveTranspose(b, x);

            // Compute residual r = b - A'x.
            AT.Multiply(-1.0, x, 1.0, r);

            Assert.IsTrue(Vector.Norm(r) < EPS);
        }

        [Test]
        public void TestSolveUnderdetermined()
        {
            // Load matrix from a file.
            var A = ResourceLoader.Get<double>("general-20x40.mat");

            Assert.IsTrue(A.RowCount < A.ColumnCount);

            // Create test data.
            var x = Helper.CreateTestVector(A.ColumnCount);
            var b = Helper.Multiply(A, x);
            var r = Vector.Clone(b);

            var qr = SparseQR.Create(A, ColumnOrdering.MinimumDegreeAtA);

            // Assuming A has full rank m, we have N(A) = n - m degrees of freedom.
            // Compute solution x with min norm(Ax - b).
            qr.Solve(b, x);

            // Compute residuals.
            A.Multiply(-1.0, x, 1.0, r);

            Assert.IsTrue(Vector.Norm(r) < EPS);
        }

        [Test]
        public void TestSolveTransposeUnderdetermined()
        {
            // Load matrix from a file.
            var A = ResourceLoader.Get<double>("general-20x40.mat");

            Assert.IsTrue(A.RowCount < A.ColumnCount);

            var AT = A.Transpose();

            // Create test data.
            var x = Helper.CreateTestVector(A.RowCount);
            var b = Helper.Multiply(AT, x);
            var r = Vector.Clone(b);

            var qr = SparseQR.Create(A, ColumnOrdering.MinimumDegreeAtA);

            qr.SolveTranspose(b, x);

            // Compute residuals.
            AT.Multiply(-1.0, x, 1.0, r);

            Assert.IsTrue(Vector.Norm(r) < EPS);
        }

        [TestCase(0, 0)]
        [TestCase(0, 5)]
        [TestCase(5, 0)]
        public void TestEmptyFactorize(int rows, int columns)
        {
            var A = new SparseMatrix(rows, columns, 0);

            var qr = SparseQR.Create(A, ColumnOrdering.MinimumDegreeAtA);

            Assert.NotNull(qr);
            Assert.IsTrue(qr.NonZerosCount == -rows);
        }
    }
}
