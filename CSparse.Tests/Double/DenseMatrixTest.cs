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
        public void TestEnumerateIndexed(int rows, int columns)
        {
            var data = MatrixHelper.LoadDense(rows, columns);

            var A = data.A;

            double sum = 0;

            A.EnumerateIndexed((i, j, a) => sum += a * a);

            Assert.That(Math.Sqrt(sum), Is.EqualTo(A.FrobeniusNorm()));
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

            var A = data.A.Clone();

            var z = Vector.Create(columns, -1.1);

            for (int i = 0; i < rows; i++)
            {
                A.SetRow(i, z);

                for (int j = 0; j < columns; j++)
                {
                    Assert.That(A.At(i, j), Is.EqualTo(-1.1));
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
                    Assert.That(A.At(i, j), Is.EqualTo(-1.1));
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

            Assert.That(actual.Values, Is.EqualTo(data.ATmB.Values).AsCollection);

            actual = A.Multiply(BT);

            Assert.That(actual.Values, Is.EqualTo(data.AmBT.Values).AsCollection);
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
            Assert.That(A.ParallelMultiply(A).Values, Is.EqualTo(A.Multiply(A).Values).AsCollection);
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
                        Assert.That(L.At(i, j), Is.EqualTo(A.At(i, j)));
                    }
                    else
                    {
                        Assert.That(L.At(i, j), Is.EqualTo(0.0));
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
                        Assert.That(L.At(i, j), Is.EqualTo(A.At(i, j)));
                    }
                    else
                    {
                        Assert.That(L.At(i, j), Is.EqualTo(0.0));
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

            Assert.That(B.RowCount, Is.EqualTo(2));
            Assert.That(B.ColumnCount, Is.EqualTo(2));

            Assert.That(B.At(0, 0), Is.EqualTo(2.0));
            Assert.That(B.At(1, 0), Is.EqualTo(2.0));
            Assert.That(B.At(0, 1), Is.EqualTo(3.0));
            Assert.That(B.At(1, 1), Is.EqualTo(3.0));

            var C = A.SubMatrix(0, 3, 0, 1);

            Assert.That(C.RowCount, Is.EqualTo(3));
            Assert.That(C.ColumnCount, Is.EqualTo(1));

            Assert.That(C.At(0, 0), Is.EqualTo(1.0));
            Assert.That(C.At(1, 0), Is.EqualTo(1.0));
            Assert.That(C.At(2, 0), Is.EqualTo(1.0));
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
                Assert.That(A.Row(i), Is.EqualTo(expected[i]).AsCollection);
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
                Assert.That(A.Row(i), Is.EqualTo(expected[i]).AsCollection);
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

            Assert.That(denseA.Equals(denseB), Is.True);
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

            Assert.That(denseA.Equals(denseB), Is.True);
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

        [TestCase(2, 2)]
        [TestCase(2, 3)]
        public void TestOfRowMajor(int rows, int columns)
        {
            var denseData = MatrixHelper.LoadDense(rows, columns);
            var denseA = denseData.A;
            var denseAT = denseA.Transpose();
            var denseB = DenseMatrix.OfRowMajor(rows, columns, denseAT.Values);

            Assert.That(denseA.Equals(denseB), Is.True);
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
                Assert.That(1.0, Is.EqualTo(A.At(i, i)));
            }
        }

        #endregion
    }
}
