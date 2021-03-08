namespace CSparse.Tests.Double
{
    using CSparse.Double;
    using NUnit.Framework;
    using System;

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
        public void TestEnumerateIndexed(int rows, int columns)
        {
            var data = MatrixHelper.LoadDense(rows, columns);

            var A = data.A;

            double sum = 0;

            A.EnumerateIndexed((i, j, a) => sum += a * a);

            Assert.AreEqual(A.FrobeniusNorm(), Math.Sqrt(sum));
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

            var A = data.A.Clone();

            var z = Vector.Create(columns, -1.1);

            for (int i = 0; i < rows; i++)
            {
                A.SetRow(i, z);

                for (int j = 0; j < columns; j++)
                {
                    Assert.AreEqual(-1.1, A.At(i, j));
                }
            }
        }

        [TestCase(2, 2)]
        [TestCase(2, 3)]
        public void TestSetColumn(int rows, int columns)
        {
            var data = MatrixHelper.LoadDense(rows, columns);

            var A = data.A.Clone();

            var z = Vector.Create(rows, -1.1);

            for (int j = 0; j < columns; j++)
            {
                A.SetColumn(j, z);

                for (int i = 0; i < rows; i++)
                {
                    Assert.AreEqual(-1.1, A.At(i, j));
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

            CollectionAssert.AreEqual(data.ATmB.Values, actual.Values);

            actual = A.Multiply(BT);

            CollectionAssert.AreEqual(data.AmBT.Values, actual.Values);
        }

        [Test]
        public void TestMatrixPointwiseMultiply()
        {
            var A = DenseMatrix.OfRowMajor(2, 3, new double[]
            {
                1.0, 2.0, 3.0,
                4.0, 5.0, 6.0
            });

            var B = DenseMatrix.OfRowMajor(2, 3, new double[]
            {
                2.0, 0.5, 2.0,
                0.5, 0.1, 0.5
            });

            var expected = DenseMatrix.OfRowMajor(2, 3, new double[]
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
            var values = new double[120 * 120];
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

        [TestCase(2, 2)]
        [TestCase(2, 3)]
        public void TestUpperTrianlge(int rows, int columns)
        {
            var data = MatrixHelper.LoadDense(rows, columns);

            var A = data.A;
            var L = A.UpperTriangle();

            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    if (i <= j)
                    {
                        Assert.AreEqual(A.At(i, j), L.At(i, j));
                    }
                    else
                    {
                        Assert.AreEqual(0.0, L.At(i, j));
                    }
                }
            }
        }

        [TestCase(2, 2)]
        [TestCase(2, 3)]
        public void TestLowerTrianlge(int rows, int columns)
        {
            var data = MatrixHelper.LoadDense(rows, columns);

            var A = data.A;
            var L = A.LowerTriangle();

            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    if (i >= j)
                    {
                        Assert.AreEqual(A.At(i, j), L.At(i, j));
                    }
                    else
                    {
                        Assert.AreEqual(0.0, L.At(i, j));
                    }
                }
            }
        }

        [Test]
        public void TestSubMatrix()
        {
            var A = DenseMatrix.OfRowMajor(3, 3, new double[]
            {
                1, 2, 3,
                1, 2, 3,
                1, 2, 3,
            });

            var B = A.SubMatrix(1, 2, 1, 2);

            Assert.AreEqual(2, B.RowCount);
            Assert.AreEqual(2, B.ColumnCount);

            Assert.AreEqual(2.0, B.At(0, 0));
            Assert.AreEqual(2.0, B.At(1, 0));
            Assert.AreEqual(3.0, B.At(0, 1));
            Assert.AreEqual(3.0, B.At(1, 1));

            var C = A.SubMatrix(0, 3, 0, 1);

            Assert.AreEqual(3, C.RowCount);
            Assert.AreEqual(1, C.ColumnCount);

            Assert.AreEqual(1.0, C.At(0, 0));
            Assert.AreEqual(1.0, C.At(1, 0));
            Assert.AreEqual(1.0, C.At(2, 0));
        }

        [Test]
        public void TestSetSubMatrix()
        {
            var A = DenseMatrix.OfRowMajor(3, 4, new double[]
            {
                0, 1, 2, 3,
                0, 1, 2, 3,
                0, 1, 2, 3
            });

            var B = DenseMatrix.OfRowMajor(2, 3, new double[]
            {
                4, 5, 6,
                7, 8, 9
            });

            A.SetSubMatrix(1, 1, B);

            var expected = new double[3][]
            {
                new double[4] { 0, 1, 2, 3 },
                new double[4] { 0, 4, 5, 6 },
                new double[4] { 0, 7, 8, 9 }
            };

            for (int i = 0; i < 3; i++)
            {
                CollectionAssert.AreEqual(expected[i], A.Row(i));
            }

            A.SetSubMatrix(0, 0, B);

            expected = new double[3][]
            {
                new double[4] { 4, 5, 6, 3 },
                new double[4] { 7, 8, 9, 6 },
                new double[4] { 0, 7, 8, 9 }
            };

            for (int i = 0; i < 3; i++)
            {
                CollectionAssert.AreEqual(expected[i], A.Row(i));
            }
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
        public void TestOfArray(int rows, int columns)
        {
            var denseData = MatrixHelper.LoadDense(rows, columns);

            var denseA = denseData.A;

            var jArray = new double[rows, columns];

            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    jArray[i, j] = denseA.At(i, j);
                }
            }

            var denseB = DenseMatrix.OfArray(jArray);

            Assert.IsTrue(denseA.Equals(denseB));
        }

        [TestCase(2, 2)]
        [TestCase(2, 3)]
        public void TestOfJaggedArray(int rows, int columns)
        {
            var denseData = MatrixHelper.LoadDense(rows, columns);

            var denseA = denseData.A;

            var jArray = new double[rows][];

            for (int i = 0; i < rows; i++)
            {
                var r = jArray[i] = new double[columns];

                for (int j = 0; j < columns; j++)
                {
                    r[j] = denseA.At(i, j);
                }
            }

            var denseB = DenseMatrix.OfJaggedArray(jArray);

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

        [TestCase(2, 2)]
        [TestCase(2, 3)]
        public void TestOfRowMajor(int rows, int columns)
        {
            var denseData = MatrixHelper.LoadDense(rows, columns);
            var denseA = denseData.A;
            var denseAT = denseA.Transpose();
            var denseB = DenseMatrix.OfRowMajor(rows, columns, denseAT.Values);

            Assert.IsTrue(denseA.Equals(denseB));
        }

        [Test]
        public void TestOfDiagonalArray()
        {
            int order = 3;

            var a = 1.0;
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
                Assert.AreEqual(A.At(i, i), 1.0);
            }
        }

        #endregion
    }
}
