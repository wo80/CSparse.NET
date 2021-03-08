namespace CSparse.Tests.Complex
{
    using CSparse.Complex;
    using NUnit.Framework;
    using System;

    using Complex = System.Numerics.Complex;

    [DefaultFloatingPointTolerance(1e-8)]
    public class DenseMatrixTest
    {
        [Test]
        public void TestConstructor()
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => new DenseMatrix(-1, 1));
            Assert.Throws<ArgumentOutOfRangeException>(() => new DenseMatrix(1, -1));
        }

        [Test]
        public void TestL1Norm()
        {
            var A = DenseMatrix.CreateDiagonal(9, 2.0);

            var actual = A.L1Norm();
            var expected = 2.0;

            Assert.AreEqual(expected, actual);

            var data = MatrixHelper.LoadDense(2, 2);

            actual = data.A.L1Norm();
            expected = 6.0;

            Assert.AreEqual(expected, actual);
        }

        [Test]
        public void TestInfinityNorm()
        {
            var A = DenseMatrix.CreateDiagonal(9, 2.0);

            var actual = A.InfinityNorm();
            var expected = 2.0;

            Assert.AreEqual(expected, actual);

            var data = MatrixHelper.LoadDense(2, 2);

            actual = data.A.InfinityNorm();
            expected = 7.0;

            Assert.AreEqual(expected, actual);
        }

        [Test]
        public void TestFrobeniusNorm()
        {
            var A = DenseMatrix.CreateDiagonal(9, 2.0);

            var actual = A.FrobeniusNorm();
            var expected = 6.0;

            Assert.AreEqual(expected, actual);

            var data = MatrixHelper.LoadDense(2, 2);

            actual = data.A.FrobeniusNorm();
            expected = 5.477225575;

            Assert.AreEqual(expected, actual);
        }

        [TestCase(2, 2)]
        [TestCase(2, 3)]
        public void TestGetRow(int rows, int columns)
        {
            var data = MatrixHelper.LoadDense(rows, columns);

            var A = data.A;

            for (int i = 0; i < rows; i++)
            {
                var y = A.Row(i);

                for (int j = 0; j < columns; j++)
                {
                    Assert.AreEqual(A.At(i, j), y[j]);
                }
            }
        }

        [TestCase(2, 2)]
        [TestCase(2, 3)]
        public void TestGetColumn(int rows, int columns)
        {
            var data = MatrixHelper.LoadDense(rows, columns);

            var A = data.A;

            for (int j = 0; j < columns; j++)
            {
                var y = A.Column(j);

                for (int i = 0; i < rows; i++)
                {
                    Assert.AreEqual(A.At(i, j), y[i]);
                }
            }
        }

        [TestCase(2, 2)]
        [TestCase(2, 3)]
        public void TestSetRow(int rows, int columns)
        {
            var data = MatrixHelper.LoadDense(rows, columns);

            var x = new Complex(-1.1, 0);

            var A = data.A.Clone();
            var z = Vector.Create(columns, x);

            for (int i = 0; i < rows; i++)
            {
                A.SetRow(i, z);

                for (int j = 0; j < columns; j++)
                {
                    Assert.AreEqual(x, A.At(i, j));
                }
            }
        }

        [TestCase(2, 2)]
        [TestCase(2, 3)]
        public void TestSetColumn(int rows, int columns)
        {
            var data = MatrixHelper.LoadDense(rows, columns);

            var x = new Complex(-1.1, 0);

            var A = data.A.Clone();
            var z = Vector.Create(rows, x);

            for (int j = 0; j < columns; j++)
            {
                A.SetColumn(j, z);

                for (int i = 0; i < rows; i++)
                {
                    Assert.AreEqual(x, A.At(i, j));
                }
            }
        }

        [TestCase(2, 2)]
        [TestCase(2, 3)]
        public void TestMatrixVectorMultiply(int rows, int columns)
        {
            var data = MatrixHelper.LoadDense(rows, columns);

            var A = data.A;
            var x = data.x;

            var actual = Vector.Create(A.RowCount, 0.0);

            A.Multiply(x, actual);

            CollectionAssert.AreEqual(data.Ax, actual);
        }

        [TestCase(2, 2)]
        [TestCase(2, 3)]
        public void TestMatrixVectorTransposeMultiply(int rows, int columns)
        {
            var data = MatrixHelper.LoadDense(rows, columns);

            var A = data.A;
            var y = data.y;

            var actual = Vector.Create(A.ColumnCount, 0.0);

            A.TransposeMultiply(y, actual);

            CollectionAssert.AreEqual(data.ATy, actual);

            var AT = data.AT;

            AT.Multiply(y, actual);

            CollectionAssert.AreEqual(data.ATy, actual);
        }

        [TestCase(2, 2)]
        [TestCase(2, 3)]
        public void TestMatrixVectorMultiplyUpdate(int rows, int columns)
        {
            var data = MatrixHelper.LoadDense(rows, columns);

            var A = data.A;
            var x = data.x;

            var actual = Vector.Create(A.RowCount, 1.0);

            A.Multiply(1.0, x, 0.0, actual);

            CollectionAssert.AreEqual(data.Ax, actual);
        }

        [TestCase(2, 2)]
        [TestCase(2, 3)]
        public void TestMatrixVectorTransposeMultiplyUpdate(int rows, int columns)
        {
            var data = MatrixHelper.LoadDense(rows, columns);

            var A = data.A;
            var y = data.y;

            var actual = Vector.Create(A.ColumnCount, 1.0);

            A.TransposeMultiply(1.0, y, 0.0, actual);

            CollectionAssert.AreEqual(data.ATy, actual);
        }

        [TestCase(2, 2)]
        [TestCase(2, 3)]
        public void TestMatrixTranspose(int rows, int columns)
        {
            var data = MatrixHelper.LoadDense(rows, columns);

            var A = data.A;
            var B = data.B;

            var actualA = A.Transpose();
            var actualB = B.Transpose();

            CollectionAssert.AreEqual(data.AT.Values, actualA.Values);
            CollectionAssert.AreEqual(data.BT.Values, actualB.Values);
        }

        [TestCase(2, 2)]
        [TestCase(2, 3)]
        public void TestMatrixSum(int rows, int columns)
        {
            var data = MatrixHelper.LoadDense(rows, columns);

            var A = data.A;
            var B = data.B;

            var actual = A.Add(B);

            CollectionAssert.AreEqual(data.ApB.Values, actual.Values);
        }

        [TestCase(2, 2)]
        [TestCase(2, 3)]
        public void TestMatrixMultiply(int rows, int columns)
        {
            var data = MatrixHelper.LoadDense(rows, columns);

            var A = data.A;
            var B = data.B;

            var AT = data.AT;
            var BT = data.BT;

            var actual = AT.Multiply(B);

            CollectionAssert.AreEqual(data.ATmB.Values, actual.Values, ComplexNumberComparer.Default);

            actual = A.Multiply(BT);

            CollectionAssert.AreEqual(data.AmBT.Values, actual.Values, ComplexNumberComparer.Default);
        }

        [Test]
        public void TestMatrixPointwiseMultiply()
        {
            var A = DenseMatrix.OfRowMajor(2, 3, new Complex[]
            {
                1.0, 2.0, 3.0,
                4.0, 5.0, 6.0
            });

            var B = DenseMatrix.OfRowMajor(2, 3, new Complex[]
            {
                2.0, 0.5, 2.0,
                0.5, 0.1, 0.5
            });

            var expected = DenseMatrix.OfRowMajor(2, 3, new Complex[]
            {
                2.0, 1.0, 6.0,
                2.0, 0.5, 3.0
            });

            var actual = new DenseMatrix(2, 3);

            A.PointwiseMultiply(B, actual);

            CollectionAssert.AreEqual(expected.Values, actual.Values);

            // Test exceptions:

            actual = new DenseMatrix(2, 2);

            Assert.Throws<ArgumentException>(() => A.PointwiseMultiply(B, actual));

            B = new DenseMatrix(2, 2);

            Assert.Throws<ArgumentException>(() => A.PointwiseMultiply(B, actual));
        }

        [Test]
        public void TestMatrixParallelMultiply()
        {
            var data = ResourceLoader.Get<double>("general-40x40.mat");
            var values = new Complex[120 * 120];
            foreach (var item in data.EnumerateIndexed())
            {
                int i = item.Item1;
                int j = item.Item2;
                for (var k = 0; k < 3; k++)
                {
                    for (var m = 0; m < 3; m++)
                    {
                        values[120 * (i + 40 * k) + j + 40 * m] = item.Item3;
                    }
                }
            }
            var A = DenseMatrix.OfColumnMajor(120, 120, values);
            CollectionAssert.AreEqual(A.Multiply(A).Values, A.ParallelMultiply(A).Values);
        }

        #region Matrix creation

        [TestCase(2, 2)]
        [TestCase(2, 3)]
        public void TestOfMatrix(int rows, int columns)
        {
            var denseData = MatrixHelper.LoadDense(rows, columns);

            var denseA = denseData.A;
            var denseB = DenseMatrix.OfMatrix(denseA);

            Assert.IsTrue(denseA.Equals(denseB));
        }

        [TestCase(2, 2)]
        [TestCase(2, 3)]
        public void TestOfIndexed(int rows, int columns)
        {
            var denseData = MatrixHelper.LoadDense(rows, columns);
            var denseA = denseData.A;
            var denseB = DenseMatrix.OfIndexed(rows, columns, denseA.EnumerateIndexed());

            Assert.IsTrue(denseA.Equals(denseB));
        }

        [TestCase(2, 2)]
        [TestCase(2, 3)]
        public void TestOfColumnMajor(int rows, int columns)
        {
            var denseData = MatrixHelper.LoadDense(rows, columns);
            var denseA = denseData.A;
            var denseB = DenseMatrix.OfColumnMajor(rows, columns, denseA.Values);

            Assert.IsTrue(denseA.Equals(denseB));
        }

        [Test]
        public void TestOfDiagonalArray()
        {
            int order = 3;

            var a = Complex.One;
            var diag = Vector.Create(order, a);

            var A = DenseMatrix.OfDiagonalArray(diag);

            for (int i = 0; i < order; i++)
            {
                Assert.AreEqual(A.At(i, i), a);
            }
        }

        [Test]
        public void TestCreateIdentity()
        {
            int order = 3;

            var A = DenseMatrix.CreateIdentity(order);

            for (int i = 0; i < order; i++)
            {
                Assert.AreEqual(A.At(i, i), Complex.One);
            }
        }

        #endregion
    }
}
