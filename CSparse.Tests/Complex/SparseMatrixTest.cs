namespace CSparse.Tests.Complex
{
    using CSparse.Complex;
    using NUnit.Framework;
    using System;

    using Complex = System.Numerics.Complex;

    public class SparseMatrixTest
    {
        [TestCase(0, 0)]
        [TestCase(0, 5)]
        [TestCase(5, 0)]
        public void TestConstructor(int rows, int columns)
        {
            var A = new SparseMatrix(rows, columns, 0);

            Assert.IsNotNull(A);
        }

        [TestCase(0, 0)]
        [TestCase(0, 5)]
        [TestCase(5, 0)]
        public void TestEmptyTranspose(int rows, int columns)
        {
            var A = new SparseMatrix(rows, columns, 0);
            var B = A.Transpose();

            Assert.IsNotNull(B);
        }

        [TestCase(0, 0)]
        [TestCase(0, 5)]
        [TestCase(5, 0)]
        public void TestEmptyAdd(int rows, int columns)
        {
            var A = new SparseMatrix(rows, columns, 0);
            var B = new SparseMatrix(rows, columns, 0);

            var C = A.Add(B);

            Assert.IsNotNull(C);
        }

        [TestCase(0, 0)]
        [TestCase(0, 5)]
        public void TestEmptyMultiply(int rows, int columns)
        {
            var A = new SparseMatrix(rows, columns, 0);
            var B = new SparseMatrix(columns, rows, 0);

            var C = A.Multiply(B);

            Assert.IsNotNull(C);
        }

        [TestCase(5, 0)]
        public void TestEmptyMultiplyInvalid(int rows, int columns)
        {
            var A = new SparseMatrix(rows, columns, 0);
            var B = new SparseMatrix(columns, rows, 0);

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
            var A = new SparseMatrix(rows, columns, 0);
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
            var A = new SparseMatrix(rows, columns, 0);
            var B = new SparseMatrix(rows, columns, 0);

            var l0 = A.InfinityNorm();
            var l1 = A.L1Norm();
            var l2 = A.FrobeniusNorm();

            Assert.IsTrue(l0 == 0.0);
            Assert.IsTrue(l1 == 0.0);
            Assert.IsTrue(l2 == 0.0);
        }

        [Test]
        public void TestMultiplyZeroMatrix()
        {
            var A = SparseMatrix.OfColumnMajor(3, 2, new Complex[6]);
            var B = SparseMatrix.OfColumnMajor(2, 4, new Complex[8]);
            var C = A.Multiply(B);

            Assert.IsTrue(C.NonZerosCount == 0);
        }
    }
}
