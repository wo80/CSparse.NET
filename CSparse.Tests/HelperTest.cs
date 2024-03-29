
namespace CSparse.Tests
{
    using CSparse.Double;
    using NUnit.Framework;

    using C = System.Numerics.Complex;

    public class HelperTest
    {
        [Test]
        public void TestZeroOfDouble()
        {
            Assert.That(Helper.ZeroOf<double>(), Is.EqualTo(0.0));
        }

        [Test]
        public void TestZeroOfComplex()
        {
            Assert.That(Helper.ZeroOf<C>(), Is.EqualTo(C.Zero));
        }

        [Test]
        public void TestTrimStorage()
        {
            var data = Double.MatrixHelper.LoadSparse(2, 2);

            var A = data.A.Clone();

            A.ColumnPointers[A.ColumnCount] = 0; // Set non-zero count to 0.

            Helper.TrimStorage(A);

            Assert.That(A.RowIndices.Length, Is.EqualTo(0));
            Assert.That(A.Values.Length, Is.EqualTo(0));
        }

        [Test]
        public void TestValidateStorage()
        {
            var ap = new int[] { 0, 3, 6 };
            var ai = new int[] { 0, 1, 2, 0, 1, 2 };
            var ax = new double[] { 0, 0, 0, 0, 0, 0 };

            var A = new SparseMatrix(3, 2, ax, ai, ap);

            Assert.That(Helper.ValidateStorage(A), Is.True);

            // Change order of column pointers.
            ap[1] = 6; ap[2] = 3;

            Assert.That(Helper.ValidateStorage(A), Is.False);

            // Revert change to column pointers.
            ap[1] = 3; ap[2] = 6;

            // Row index larger than number of rows.
            ai[2] = 3;

            Assert.That(Helper.ValidateStorage(A), Is.False);

            // Change order of row indices.
            ai[1] = 2; ai[2] = 1;

            Assert.That(Helper.ValidateStorage(A), Is.True);
            Assert.That(Helper.ValidateStorage(A, true), Is.False);
        }
    }
}
