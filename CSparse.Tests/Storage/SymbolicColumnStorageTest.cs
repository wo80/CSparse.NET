
namespace CSparse.Tests.Storage
{
    using CSparse.Storage;
    using NUnit.Framework;
    using System;

    public class SymbolicColumnStorageTest
    {
        [TestCase(0, 0)]
        [TestCase(0, 5)]
        [TestCase(5, 0)]
        public void TestConstructor(int rows, int columns)
        {
            var A = new SymbolicColumnStorage(rows, columns, 0, true);

            Assert.IsNotNull(A);
        }

        [TestCase(0, 0)]
        [TestCase(0, 5)]
        [TestCase(5, 0)]
        public void TestEmptyTranspose(int rows, int columns)
        {
            var A = new SymbolicColumnStorage(rows, columns, 0, true);
            var B = A.Transpose();

            Assert.IsNotNull(B);
        }

        [TestCase(0, 0)]
        [TestCase(0, 5)]
        [TestCase(5, 0)]
        public void TestEmptyAdd(int rows, int columns)
        {
            var A = new SymbolicColumnStorage(rows, columns, 0, true);
            var B = new SymbolicColumnStorage(rows, columns, 0, true);

            var C = A.Add(B);

            Assert.IsNotNull(C);
        }
    }
}
