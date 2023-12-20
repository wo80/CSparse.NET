
namespace CSparse.Tests
{
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
    }
}
