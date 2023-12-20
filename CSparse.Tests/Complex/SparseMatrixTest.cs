namespace CSparse.Tests.Complex
{
    using CSparse.Complex;
    using CSparse.Storage;
    using NUnit.Framework;
    using System;

    using Complex = System.Numerics.Complex;

    [DefaultFloatingPointTolerance(1e-8)]
    public class SparseMatrixTest
    {
        #region Test empty matrix

        [TestCase(0, 0)]
        [TestCase(0, 5)]
        [TestCase(5, 0)]
        public void TestConstructor(int rows, int columns)
        {
            var A = new SparseMatrix(rows, columns, 0);

            Assert.That(A, Is.Not.Null);
        }

        [TestCase(0, 0)]
        [TestCase(0, 5)]
        [TestCase(5, 0)]
        public void TestEmptyTranspose(int rows, int columns)
        {
            var A = new SparseMatrix(rows, columns, 0);
            var B = A.Transpose();

            Assert.That(B, Is.Not.Null);
        }

        [TestCase(0, 0)]
        [TestCase(0, 5)]
        [TestCase(5, 0)]
        public void TestEmptyAdd(int rows, int columns)
        {
            var A = new SparseMatrix(rows, columns, 0);
            var B = new SparseMatrix(rows, columns, 0);

            var C = A.Add(B);

            Assert.That(C, Is.Not.Null);
        }

        [TestCase(0, 0)]
        [TestCase(0, 5)]
        public void TestEmptyMultiply(int rows, int columns)
        {
            var A = new SparseMatrix(rows, columns, 0);
            var B = new SparseMatrix(columns, rows, 0);

            var C = A.Multiply(B);

            Assert.That(C, Is.Not.Null);
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

            Assert.That(Vector.Norm(y.Length, y), Is.EqualTo(0.0));
        }

        [TestCase(0, 0)]
        [TestCase(0, 5)]
        [TestCase(5, 0)]
        public void TestEmptyNorm(int rows, int columns)
        {
            var A = new SparseMatrix(rows, columns, 0);

            var l0 = A.InfinityNorm();
            var l1 = A.L1Norm();
            var l2 = A.FrobeniusNorm();

            Assert.That(l0, Is.EqualTo(0.0));
            Assert.That(l1, Is.EqualTo(0.0));
            Assert.That(l2, Is.EqualTo(0.0));
        }

        [Test]
        public void TestMultiplyZeroMatrix()
        {
            var A = SparseMatrix.OfColumnMajor(3, 2, new Complex[6]);
            var B = SparseMatrix.OfColumnMajor(2, 4, new Complex[8]);
            var C = A.Multiply(B);

            Assert.That(C.NonZerosCount, Is.EqualTo(0));
        }

        #endregion

        [Test]
        public void TestConstructor()
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => new SparseMatrix(-1, 1));
            Assert.Throws<ArgumentOutOfRangeException>(() => new SparseMatrix(1, -1));
            Assert.Throws<ArgumentOutOfRangeException>(() => new SparseMatrix(1, 1, -1));
        }

        [Test]
        public void TestDropZeros()
        {
            var data = MatrixHelper.LoadSparse(2, 2);

            var A = data.A.Clone();

            A.Values[0] = Complex.Zero;

            int nnz = A.DropZeros();

            Assert.That(3, Is.EqualTo(nnz));
            Assert.That(3, Is.EqualTo(A.NonZerosCount));
        }

        [Test]
        public void TestKeep()
        {
            var data = MatrixHelper.LoadSparse(2, 2);

            var A = data.A.Clone();

            // Keep strict upper part of the matrix.
            int nnz = A.Keep((i, j, a) => i < j);

            Assert.That(1, Is.EqualTo(nnz));
            Assert.That(1, Is.EqualTo(A.NonZerosCount));
        }

        [Test]
        public void TestL1Norm()
        {
            var A = SparseMatrix.CreateDiagonal(9, 2.0);

            var actual = A.L1Norm();
            var expected = 2.0;

            Assert.That(actual, Is.EqualTo(expected));

            var data = MatrixHelper.LoadSparse(2, 2);

            actual = data.A.L1Norm();
            expected = 6.0;

            Assert.That(actual, Is.EqualTo(expected));
        }

        [Test]
        public void TestInfinityNorm()
        {
            var A = SparseMatrix.CreateDiagonal(9, 2.0);

            var actual = A.InfinityNorm();
            var expected = 2.0;

            Assert.That(actual, Is.EqualTo(expected));

            var data = MatrixHelper.LoadSparse(2, 2);

            actual = data.A.InfinityNorm();
            expected = 7.0;

            Assert.That(actual, Is.EqualTo(expected));
        }

        [Test]
        public void TestFrobeniusNorm()
        {
            var A = SparseMatrix.CreateDiagonal(9, 2.0);

            var actual = A.FrobeniusNorm();
            var expected = 6.0;

            Assert.That(actual, Is.EqualTo(expected));

            var data = MatrixHelper.LoadSparse(2, 2);

            actual = data.A.FrobeniusNorm();
            expected = 5.477225575;

            Assert.That(actual, Is.EqualTo(expected));
        }

        [TestCase(2, 2)]
        [TestCase(2, 3)]
        public void TestEnumerateIndexed(int rows, int columns)
        {
            var data = MatrixHelper.LoadSparse(rows, columns);

            var A = data.A;

            double sum = 0;

            A.EnumerateIndexed((i, j, a) => sum += a.Real * a.Real + a.Imaginary * a.Imaginary);

            Assert.That(Math.Sqrt(sum), Is.EqualTo(A.FrobeniusNorm()));
        }

        [TestCase(2, 2)]
        [TestCase(2, 3)]
        public void TestGetRow(int rows, int columns)
        {
            var data = MatrixHelper.LoadSparse(rows, columns);

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
            var data = MatrixHelper.LoadSparse(rows, columns);

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
        public void TestMatrixVectorMultiply(int rows, int columns)
        {
            var data = MatrixHelper.LoadSparse(rows, columns);

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
            var data = MatrixHelper.LoadSparse(rows, columns);

            var A = data.A;
            var y = data.y;

            var actual = Vector.Create(A.ColumnCount, 0.0);

            A.TransposeMultiply(y, actual);

            Assert.That(actual, Is.EqualTo(data.ATy).AsCollection);

            Vector.Clear(actual);

            var AT = data.AT;

            AT.Multiply(y, actual);

            Assert.That(actual, Is.EqualTo(data.ATy).AsCollection);
        }

        [TestCase(2, 2)]
        [TestCase(2, 3)]
        public void TestMatrixVectorMultiplyUpdate(int rows, int columns)
        {
            var data = MatrixHelper.LoadSparse(rows, columns);

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
            var data = MatrixHelper.LoadSparse(rows, columns);

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
            var data = MatrixHelper.LoadSparse(rows, columns);

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
            var data = MatrixHelper.LoadSparse(rows, columns);

            var A = data.A;
            var B = data.B;

            var actual = A.Add(B);

            Assert.That(actual.Values, Is.EqualTo(data.ApB.Values).AsCollection);
        }

        [TestCase(2, 2)]
        [TestCase(2, 3)]
        public void TestMatrixMultiply(int rows, int columns)
        {
            var data = MatrixHelper.LoadSparse(rows, columns);

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
        public void TestMatrixParallelMultiply()
        {
            var data = ResourceLoader.Get<double>("general-40x40.mat");
            var acs = new CoordinateStorage<double>(40, 800, 40 * 800);
            var bcs = new CoordinateStorage<double>(800, 40, 800 * 40);
            // This just exceeds min_total_ops in ParallelMultiply
            foreach (var item in data.EnumerateIndexed())
            {
                int i = item.Item1;
                int j = item.Item2;
                for (var k = 0; k < 20; k++)
                {
                    acs.At(i, j + 40 * k, item.Item3);
                    bcs.At(i + 40 * k, j, item.Item3);
                }
            }
            var A = CompressedColumnStorage<double>.OfIndexed(acs);
            var B = CompressedColumnStorage<double>.OfIndexed(bcs);

            var sResult = A.Multiply(B);
            var pResult = A.ParallelMultiply(B);

            Assert.That(pResult.Values, Is.EqualTo(sResult.Values).AsCollection);
        }

        [TestCase(2, 2)]
        [TestCase(2, 3)]
        public void TestMatrixPermuteColumns(int rows, int columns)
        {
            var p = Permutation.Create(columns, -1);

            var data = MatrixHelper.LoadSparse(rows, columns);

            var A = data.A;
            var Ap = A.PermuteColumns(p);

            var actualColumn = new Complex[rows];
            var expectedColumn = new Complex[rows];

            for (int i = 0; i < columns; i++)
            {
                A.Column(p[i], expectedColumn);
                Ap.Column(i, actualColumn);

                Assert.That(actualColumn, Is.EqualTo(expectedColumn).AsCollection);
            }
        }

        [TestCase(2, 2)]
        [TestCase(2, 3)]
        public void TestMatrixPermuteRows(int rows, int columns)
        {
            var p = Permutation.Create(rows, -1);

            var data = MatrixHelper.LoadSparse(rows, columns);

            var A = data.A;
            var Ap = A.Clone();

            var actualRow = new Complex[columns];
            var expectedRow = new Complex[columns];

            Ap.PermuteRows(p);

            for (int i = 0; i < rows; i++)
            {
                A.Row(p[i], expectedRow);
                Ap.Row(i, actualRow);

                Assert.That(actualRow, Is.EqualTo(expectedRow).AsCollection);
            }
        }

        #region Matrix creation

        [TestCase(2, 2)]
        [TestCase(2, 3)]
        public void TestOfMatrix(int rows, int columns)
        {
            var sparseData = MatrixHelper.LoadSparse(rows, columns);

            var sparseA = sparseData.A;
            var sparseB = SparseMatrix.OfMatrix(sparseA);

            Assert.That(sparseA.Equals(sparseB), Is.True);

            var denseData = MatrixHelper.LoadDense(rows, columns);

            sparseB = SparseMatrix.OfMatrix(denseData.A);

            Assert.That(sparseA.Equals(sparseB), Is.True);
        }

        [TestCase(2, 2)]
        [TestCase(2, 3)]
        public void TestOfIndexed(int rows, int columns)
        {
            var sparseData = MatrixHelper.LoadSparse(rows, columns);

            var sparseA = sparseData.A;
            var sparseB = SparseMatrix.OfIndexed(rows, columns, sparseA.EnumerateIndexed());

            Assert.That(sparseA.Equals(sparseB), Is.True);
        }

        [TestCase(2, 2)]
        [TestCase(2, 3)]
        public void TestOfIndexed_Coo(int rows, int columns)
        {
            var sparseData = MatrixHelper.LoadSparse(rows, columns);

            var sparseA = sparseData.A;

            var coo = new CoordinateStorage<Complex>(rows, columns, sparseA.NonZerosCount);

            foreach (var a in sparseA.EnumerateIndexed())
            {
                coo.At(a.Item1, a.Item2, a.Item3);
            }

            var sparseB = SparseMatrix.OfIndexed(coo);

            Assert.That(sparseA.Equals(sparseB), Is.True);

            var sparseC = SparseMatrix.OfIndexed(coo, true);

            Assert.That(sparseA.Equals(sparseC), Is.True);
            Assert.That(coo.Values, Is.Null);
            Assert.That(coo.RowIndices, Is.Null);
            Assert.That(coo.ColumnIndices, Is.Null);
        }

        [TestCase(2, 2)]
        [TestCase(2, 3)]
        public void TestOfColumnMajor(int rows, int columns)
        {
            var denseData = MatrixHelper.LoadDense(rows, columns);
            var denseA = denseData.A;

            var sparseData = MatrixHelper.LoadSparse(rows, columns);
            var sparseA = sparseData.A;

            var sparseB = SparseMatrix.OfColumnMajor(rows, columns, denseA.Values);

            Assert.That(sparseA.Equals(sparseB), Is.True);
        }

        [TestCase(2, 2)]
        [TestCase(2, 3)]
        public void TestOfArray(int rows, int columns)
        {
            var sparseData = MatrixHelper.LoadSparse(rows, columns);

            var sparseA = sparseData.A;

            var array = new Complex[rows, columns];

            foreach (var a in sparseA.EnumerateIndexed())
            {
                array[a.Item1, a.Item2] = a.Item3;
            }

            var sparseB = SparseMatrix.OfArray(array);

            Assert.That(sparseA.Equals(sparseB), Is.True);
        }

        [TestCase(2, 2)]
        [TestCase(2, 3)]
        public void TestOfJaggedArray(int rows, int columns)
        {
            var sparseData = MatrixHelper.LoadSparse(rows, columns);

            var sparseA = sparseData.A;

            var array = new Complex[rows][];

            for (int i = 0; i < rows; i++)
            {
                var values = new Complex[columns];

                for (int j = 0; j < columns; j++)
                {
                    values[j] = sparseA.At(i, j);
                }

                array[i] = values;
            }

            var sparseB = SparseMatrix.OfJaggedArray(array);

            Assert.That(sparseA.Equals(sparseB), Is.True);
        }

        [Test]
        public void TestOfDiagonalArray()
        {
            int order = 3;

            var a = Complex.One;
            var diag = Vector.Create(order, a);

            var A = SparseMatrix.OfDiagonalArray(diag);

            for (int i = 0; i < order; i++)
            {
                Assert.That(a, Is.EqualTo(A.At(i, i)));
            }
        }

        [Test]
        public void TestCreateIdentity()
        {
            int order = 3;

            var A = SparseMatrix.CreateIdentity(order);

            for (int i = 0; i < order; i++)
            {
                Assert.That(Complex.One, Is.EqualTo(A.At(i, i)));
            }
        }

        #endregion
    }
}
