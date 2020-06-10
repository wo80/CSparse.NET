
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
    }
}
