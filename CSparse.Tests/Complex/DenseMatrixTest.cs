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

            Assert.That(actual, Is.EqualTo(expected));

            var data = MatrixHelper.LoadDense(2, 2);

            actual = data.A.L1Norm();
            expected = 6.0;

            Assert.That(actual, Is.EqualTo(expected));
        }

        [Test]
        public void TestInfinityNorm()
        {
            var A = DenseMatrix.CreateDiagonal(9, 2.0);

            var actual = A.InfinityNorm();
            var expected = 2.0;

            Assert.That(actual, Is.EqualTo(expected));

            var data = MatrixHelper.LoadDense(2, 2);

            actual = data.A.InfinityNorm();
            expected = 7.0;

            Assert.That(actual, Is.EqualTo(expected));
        }

        [Test]
        public void TestFrobeniusNorm()
        {
            var A = DenseMatrix.CreateDiagonal(9, 2.0);

            var actual = A.FrobeniusNorm();
            var expected = 6.0;

            Assert.That(actual, Is.EqualTo(expected));

            var data = MatrixHelper.LoadDense(2, 2);

            actual = data.A.FrobeniusNorm();
            expected = 5.477225575;

            Assert.That(actual, Is.EqualTo(expected));
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
                    Assert.That(y[j], Is.EqualTo(A.At(i, j)));
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
                    Assert.That(y[i], Is.EqualTo(A.At(i, j)));
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
                    Assert.That(A.At(i, j), Is.EqualTo(x));
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
                    Assert.That(A.At(i, j), Is.EqualTo(x));
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

            Assert.That(actual, Is.EqualTo(data.Ax).AsCollection);
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

            Assert.That(actual, Is.EqualTo(data.ATy).AsCollection);

            var AT = data.AT;

            AT.Multiply(y, actual);

            Assert.That(actual, Is.EqualTo(data.ATy).AsCollection);
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

            Assert.That(actual, Is.EqualTo(data.Ax).AsCollection);
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

            Assert.That(actual, Is.EqualTo(data.ATy).AsCollection);
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

            Assert.That(actualA.Values, Is.EqualTo(data.AT.Values).AsCollection);
            Assert.That(actualB.Values, Is.EqualTo(data.BT.Values).AsCollection);
        }

        [TestCase(2, 2)]
        [TestCase(2, 3)]
        public void TestMatrixSum(int rows, int columns)
        {
            var data = MatrixHelper.LoadDense(rows, columns);

            var A = data.A;
            var B = data.B;

            var actual = A.Add(B);

            Assert.That(actual.Values, Is.EqualTo(data.ApB.Values).AsCollection);
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

            Assert.That(actual.Values, Is.EqualTo(data.ATmB.Values).Using<Complex>(ComplexNumberComparer.Default));

            actual = A.Multiply(BT);

            Assert.That(actual.Values, Is.EqualTo(data.AmBT.Values).Using<Complex>(ComplexNumberComparer.Default));
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

            Assert.That(actual.Values, Is.EqualTo(expected.Values).AsCollection);

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
            Assert.That(A.ParallelMultiply(A).Values, Is.EqualTo(A.Multiply(A).Values).AsCollection);
        }

        #region Matrix creation

        [TestCase(2, 2)]
        [TestCase(2, 3)]
        public void TestOfMatrix(int rows, int columns)
        {
            var denseData = MatrixHelper.LoadDense(rows, columns);

            var denseA = denseData.A;
            var denseB = DenseMatrix.OfMatrix(denseA);

            Assert.That(denseA.Equals(denseB), Is.True);
        }

        [TestCase(2, 2)]
        [TestCase(2, 3)]
        public void TestOfIndexed(int rows, int columns)
        {
            var denseData = MatrixHelper.LoadDense(rows, columns);
            var denseA = denseData.A;
            var denseB = DenseMatrix.OfIndexed(rows, columns, denseA.EnumerateIndexed());

            Assert.That(denseA.Equals(denseB), Is.True);
        }

        [TestCase(2, 2)]
        [TestCase(2, 3)]
        public void TestOfColumnMajor(int rows, int columns)
        {
            var denseData = MatrixHelper.LoadDense(rows, columns);
            var denseA = denseData.A;
            var denseB = DenseMatrix.OfColumnMajor(rows, columns, denseA.Values);

            Assert.That(denseA.Equals(denseB), Is.True);
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
                Assert.That(a, Is.EqualTo(A.At(i, i)));
            }
        }

        [Test]
        public void TestCreateIdentity()
        {
            int order = 3;

            var A = DenseMatrix.CreateIdentity(order);

            for (int i = 0; i < order; i++)
            {
                Assert.That(Complex.One, Is.EqualTo(A.At(i, i)));
            }
        }

        #endregion
    }
}
