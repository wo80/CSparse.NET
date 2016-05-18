
namespace CSparse.Tests.Double
{
    using CSparse.Double;
    using NUnit.Framework;
    using System;

    public class CompressedColumnStorageTest
    {
        [TestCase(0, 0)]
        [TestCase(0, 5)]
        [TestCase(5, 0)]
        public void TestConstructor(int rows, int columns)
        {
            var A = MatrixHelper.Load(rows, columns);

            Assert.IsNotNull(A);
        }

        [TestCase(0, 0)]
        [TestCase(0, 5)]
        [TestCase(5, 0)]
        public void TestEmptyTranspose(int rows, int columns)
        {
            var A = MatrixHelper.Load(rows, columns);
            var B = A.Transpose();

            Assert.IsNotNull(B);
        }

        [TestCase(0, 0)]
        [TestCase(0, 5)]
        [TestCase(5, 0)]
        public void TestEmptyAdd(int rows, int columns)
        {
            var A = MatrixHelper.Load(rows, columns);
            var B = MatrixHelper.Load(rows, columns);

            var C = A.Add(B);

            Assert.IsNotNull(C);
        }

        [TestCase(0, 0)]
        [TestCase(0, 5)]
        public void TestEmptyMultiply(int rows, int columns)
        {
            var A = MatrixHelper.Load(rows, columns);
            var B = MatrixHelper.Load(columns, rows);

            var C = A.Multiply(B);

            Assert.IsNotNull(C);
        }

        [TestCase(5, 0)]
        public void TestEmptyMultiplyInvalid(int rows, int columns)
        {
            var A = MatrixHelper.Load(rows, columns);
            var B = MatrixHelper.Load(columns, rows);

            Assert.Throws<Exception>(() =>
            {
                var C = A.Multiply(B);
            });
        }

        [TestCase(0, 0)]
        [TestCase(0, 5)]
        [TestCase(5, 0)]
        public void TestEmptyVectorMultiply(int rows, int columns)
        {
            var A = MatrixHelper.Load(rows, columns);
            var x = Vector.Create(columns, 1.0);
            var y = Vector.Create(rows, 0.0);

            A.Multiply(x, y);

            Assert.IsTrue(Vector.Norm(y) == 0.0);
        }

        [TestCase(0, 0)]
        [TestCase(0, 5)]
        [TestCase(5, 0)]
        public void TestEmptyNorm(int rows, int columns)
        {
            var A = MatrixHelper.Load(rows, columns);
            var B = MatrixHelper.Load(rows, columns);

            var l0 = A.Norm(0);
            var l1 = A.Norm(1);
            var l2 = A.Norm(2);

            Assert.IsTrue(l0 == 0.0);
            Assert.IsTrue(l1 == 0.0);
            Assert.IsTrue(l2 == 0.0);
        }
    }
}
