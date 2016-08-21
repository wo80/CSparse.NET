
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

        [Test]
        public void TestJaggedArrayToCSC(int n)
        {
            double[][] array = new double[3][];
            for (int i = 0; i < 3; i++)
                array[i] = new double[3];

            array[0][1] = 2.0; array[1][0] = 2.0;
            array[0][2] = 4.0; array[2][0] = 4.0;
            array[2][2] = 6.0;

            var csc = CSparse.Converter.ToCompressedColumnStorage(array, false);
            Assert.IsTrue(csc.IsSymmetric());

            double[] x = new[] { 1.0, 2.0, 3.0 };
            double[] y = new double[3];
            csc.Multiply(x, y);

            Assert.AreEqual(y[0], 16.0 ,1e-6);
            Assert.AreEqual(y[1], 2.0 , 1e-6);
            Assert.AreEqual(y[2], 22.0, 1e-6);
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
