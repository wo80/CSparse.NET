
namespace CSparse.Tests
{
    using NUnit.Framework;

    using C = System.Numerics.Complex;

    public class HelperTest
    {
        [Test]
        public void TestZeroOfDouble()
        {
            Assert.AreEqual(0.0, Helper.ZeroOf<double>());
        }

        [Test]
        public void TestZeroOfComplex()
        {
            Assert.AreEqual(C.Zero, Helper.ZeroOf<C>());
        }

        [Test]
        public void TestTrimStorage()
        {
            var data = Double.MatrixHelper.LoadSparse(2, 2);

            var A = data.A;

            A.ColumnPointers[A.ColumnCount] = 0; // Set non-zero count to 0.

            Helper.TrimStorage(A);

            Assert.AreEqual(0, A.RowIndices.Length);
            Assert.AreEqual(0, A.Values.Length);
        }
    }
}
