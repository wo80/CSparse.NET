namespace CSparse.Tests.Complex.Factorization
{
    using CSparse.Complex;
    using CSparse.Complex.Factorization;
    using NUnit.Framework;
    using Complex = System.Numerics.Complex;

    public class SparseQRTest
    {
        private const double EPS = 1.0e-6;

        [Test]
        public void TestSolve()
        {
            // Load matrix from a file.
            var A = ResourceLoader.Get<Complex>("general-40x40.mat");

            Assert.That(A.ColumnCount, Is.EqualTo(A.RowCount));

            // Create test data.
            var x = Helper.CreateTestVector(A.ColumnCount);
            var b = Helper.Multiply(A, x);
            var r = Vector.Clone(b);

            var qr = SparseQR.Create(A, ColumnOrdering.MinimumDegreeAtA);

            // Solve Ax = b.
            qr.Solve(b, x);

            // Compute residual r = b - Ax.
            A.Multiply(-1.0, x, 1.0, r);

            Assert.That(Vector.Norm(r.Length, r) < EPS, Is.True);
        }

        [Test]
        public void TestSolveTranspose()
        {
            // Load matrix from a file.
            var A = ResourceLoader.Get<Complex>("general-40x40.mat");

            Assert.That(A.ColumnCount, Is.EqualTo(A.RowCount));

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

            Assert.That(Vector.Norm(r.Length, r) < EPS, Is.True);
        }

        [Test]
        public void TestSolveOverdetermined()
        {
            // Load matrix from a file.
            var A = ResourceLoader.Get<Complex>("general-40x20.mat");

            Assert.That(A.RowCount > A.ColumnCount, Is.True);

            // Create test data.
            var x = Helper.CreateTestVector(A.ColumnCount);
            var b = Helper.Multiply(A, x);
            var r = Vector.Clone(b);

            var qr = SparseQR.Create(A, ColumnOrdering.MinimumDegreeAtA);

            // Compute min norm(Ax - b).
            qr.Solve(b, x);

            // Compute residual r = b - Ax.
            A.Multiply(-1.0, x, 1.0, r);

            Assert.That(Vector.Norm(r.Length, r) < EPS, Is.True);
        }

        [Test]
        public void TestSolveTransposeOverdetermined()
        {
            // Load matrix from a file.
            var A = ResourceLoader.Get<Complex>("general-40x20.mat");

            Assert.That(A.RowCount > A.ColumnCount, Is.True);

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

            Assert.That(Vector.Norm(r.Length, r) < EPS, Is.True);
        }

        [Test]
        public void TestSolveUnderdetermined()
        {
            // Load matrix from a file.
            var A = ResourceLoader.Get<Complex>("general-20x40.mat");

            Assert.That(A.RowCount < A.ColumnCount, Is.True);

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

            Assert.That(Vector.Norm(r.Length, r) < EPS, Is.True);
        }

        [Test]
        public void TestSolveTransposeUnderdetermined()
        {
            // Load matrix from a file.
            var A = ResourceLoader.Get<Complex>("general-20x40.mat");

            Assert.That(A.RowCount < A.ColumnCount, Is.True);

            var AT = A.Transpose();

            // Create test data.
            var x = Helper.CreateTestVector(A.RowCount);
            var b = Helper.Multiply(AT, x);
            var r = Vector.Clone(b);

            var qr = SparseQR.Create(A, ColumnOrdering.MinimumDegreeAtA);

            qr.SolveTranspose(b, x);

            // Compute residuals.
            AT.Multiply(-1.0, x, 1.0, r);

            Assert.That(Vector.Norm(r.Length, r) < EPS, Is.True);
        }

        [TestCase(0, 0)]
        [TestCase(0, 5)]
        [TestCase(5, 0)]
        public void TestEmptyFactorize(int rows, int columns)
        {
            var A = new SparseMatrix(rows, columns, 0);

            var qr = SparseQR.Create(A, ColumnOrdering.MinimumDegreeAtA);

            Assert.That(qr, Is.Not.Null);
            Assert.That(qr.NonZerosCount == -rows, Is.True);
        }
    }
}
