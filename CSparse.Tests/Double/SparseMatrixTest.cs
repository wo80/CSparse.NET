namespace CSparse.Tests.Double
{
    using CSparse.Double;
    using CSparse.Storage;
    using NUnit.Framework;
    using System;

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

            Assert.IsTrue(Vector.Norm(y.Length, y) == 0.0);
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

            Assert.IsTrue(l0 == 0.0);
            Assert.IsTrue(l1 == 0.0);
            Assert.IsTrue(l2 == 0.0);
        }

        [Test]
        public void TestMultiplyZeroMatrix()
        {
            var A = SparseMatrix.OfColumnMajor(3, 2, new double[6]);
            var B = SparseMatrix.OfColumnMajor(2, 4, new double[8]);
            var C = A.Multiply(B);

            Assert.IsTrue(C.NonZerosCount == 0);
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

            A.Values[0] = 0.0;

            int nnz = A.DropZeros();

            Assert.AreEqual(nnz, 3);
            Assert.AreEqual(A.NonZerosCount, 3);
        }

        [Test]
        public void TestKeep()
        {
            var data = MatrixHelper.LoadSparse(2, 2);

            var A = data.A.Clone();

            // Keep strict upper part of the matrix.
            int nnz = A.Keep((i, j, a) => i < j);

            Assert.AreEqual(nnz, 1);
            Assert.AreEqual(A.NonZerosCount, 1);
        }

        [Test]
        public void TestL1Norm()
        {
            var A = SparseMatrix.CreateDiagonal(9, 2.0);

            var actual = A.L1Norm();
            var expected = 2.0;

            Assert.AreEqual(expected, actual);

            var data = MatrixHelper.LoadSparse(2, 2);

            actual = data.A.L1Norm();
            expected = 6.0;

            Assert.AreEqual(expected, actual);
        }

        [Test]
        public void TestInfinityNorm()
        {
            var A = SparseMatrix.CreateDiagonal(9, 2.0);

            var actual = A.InfinityNorm();
            var expected = 2.0;

            Assert.AreEqual(expected, actual);

            var data = MatrixHelper.LoadSparse(2, 2);

            actual = data.A.InfinityNorm();
            expected = 7.0;

            Assert.AreEqual(expected, actual);
        }

        [Test]
        public void TestFrobeniusNorm()
        {
            var A = SparseMatrix.CreateDiagonal(9, 2.0);

            var actual = A.FrobeniusNorm();
            var expected = 6.0;

            Assert.AreEqual(expected, actual);

            var data = MatrixHelper.LoadSparse(2, 2);

            actual = data.A.FrobeniusNorm();
            expected = 5.477225575;

            Assert.AreEqual(expected, actual);
        }

        [TestCase(2, 2)]
        [TestCase(2, 3)]
        public void TestEnumerateIndexed(int rows, int columns)
        {
            var data = MatrixHelper.LoadSparse(rows, columns);

            var A = data.A;

            double sum = 0;

            A.EnumerateIndexed((i, j, a) => sum += a * a);

            Assert.AreEqual(A.FrobeniusNorm(), Math.Sqrt(sum));
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
                    Assert.AreEqual(A.At(i, j), y[j]);
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
                    Assert.AreEqual(A.At(i, j), y[i]);
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

            CollectionAssert.AreEqual(data.Ax, actual);
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

            CollectionAssert.AreEqual(data.ATy, actual);

            Vector.Clear(actual);

            var AT = data.AT;

            AT.Multiply(y, actual);

            CollectionAssert.AreEqual(data.ATy, actual);
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

            CollectionAssert.AreEqual(data.Ax, actual);
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

            CollectionAssert.AreEqual(data.ATy, actual);
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

            CollectionAssert.AreEqual(data.AT.Values, actualA.Values);
            CollectionAssert.AreEqual(data.BT.Values, actualB.Values);
        }

        [TestCase(2, 2)]
        [TestCase(2, 3)]
        public void TestMatrixSum(int rows, int columns)
        {
            var data = MatrixHelper.LoadSparse(rows, columns);

            var A = data.A;
            var B = data.B;

            var actual = A.Add(B);

            CollectionAssert.AreEqual(data.ApB.Values, actual.Values);
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

            CollectionAssert.AreEqual(data.ATmB.Values, actual.Values);

            actual = A.Multiply(BT);

            CollectionAssert.AreEqual(data.AmBT.Values, actual.Values);
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
            CollectionAssert.AreEqual(A.Multiply(B).Values, A.ParallelMultiply(B).Values);
        }

        [TestCase(2, 2)]
        [TestCase(2, 3)]
        public void TestMatrixPermuteColumns(int rows, int columns)
        {
            var p = Permutation.Create(columns, -1);

            var data = MatrixHelper.LoadSparse(rows, columns);

            var A = data.A;
            var Ap = A.PermuteColumns(p);

            var actualColumn = new double[rows];
            var expectedColumn = new double[rows];

            for (int i = 0; i < columns; i++)
            {
                A.Column(p[i], expectedColumn);
                Ap.Column(i, actualColumn);

                CollectionAssert.AreEqual(expectedColumn, actualColumn);
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

            var actualRow = new double[columns];
            var expectedRow = new double[columns];

            Ap.PermuteRows(p);

            for (int i = 0; i < rows; i++)
            {
                A.Row(p[i], expectedRow);
                Ap.Row(i, actualRow);

                CollectionAssert.AreEqual(expectedRow, actualRow);
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

            Assert.IsTrue(sparseA.Equals(sparseB));

            var denseData = MatrixHelper.LoadDense(rows, columns);

            sparseB = SparseMatrix.OfMatrix(denseData.A);

            Assert.IsTrue(sparseA.Equals(sparseB));
        }

        [TestCase(2, 2)]
        [TestCase(2, 3)]
        public void TestOfIndexed(int rows, int columns)
        {
            var sparseData = MatrixHelper.LoadSparse(rows, columns);

            var sparseA = sparseData.A;
            var sparseB = SparseMatrix.OfIndexed(rows, columns, sparseA.EnumerateIndexed());

            Assert.IsTrue(sparseA.Equals(sparseB));
        }

        [TestCase(2, 2)]
        [TestCase(2, 3)]
        public void TestOfIndexed_Coo(int rows, int columns)
        {
            var sparseData = MatrixHelper.LoadSparse(rows, columns);

            var sparseA = sparseData.A;

            var coo = new CoordinateStorage<double>(rows, columns, sparseA.NonZerosCount);

            foreach (var a in sparseA.EnumerateIndexed())
            {
                coo.At(a.Item1, a.Item2, a.Item3);
            }

            var sparseB = SparseMatrix.OfIndexed(coo);

            Assert.IsTrue(sparseA.Equals(sparseB));

            var sparseC = SparseMatrix.OfIndexed(coo, true);

            Assert.IsTrue(sparseA.Equals(sparseC));
            Assert.IsNull(coo.Values);
            Assert.IsNull(coo.RowIndices);
            Assert.IsNull(coo.ColumnIndices);
        }

        [Test]
        public void TestOfIndexed_Coo_InPlace()
        {
            // rows > columns

            var coo = new CoordinateStorage<double>(10, 5, 3);

            coo.At(0, 0, 1.0);
            coo.At(1, 1, 1.0);
            coo.At(4, 4, 1.0);

            var A = SparseMatrix.OfIndexed(coo, true);

            Assert.AreEqual(3, A.NonZerosCount);

            // rows < columns

            coo = new CoordinateStorage<double>(5, 10, 3);

            coo.At(4, 4, 1.0);
            coo.At(1, 1, 1.0);
            coo.At(0, 0, 1.0);
            coo.At(4, 4, -1.0); // Results in one explicit zero entry.

            A = SparseMatrix.OfIndexed(coo, true);

            Assert.AreEqual(3, A.NonZerosCount);

            A.DropZeros();

            Assert.AreEqual(2, A.NonZerosCount);
        }

        [Test]
        public void TestOfIndexed_Empty()
        {
            var coord = new CoordinateStorage<double>(0, 0, 0);

            var sparseA = new SparseMatrix(0, 0, 0);
            var sparseB = SparseMatrix.OfIndexed(coord);

            Assert.IsTrue(sparseA.Equals(sparseB));
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

            Assert.IsTrue(sparseA.Equals(sparseB));
        }

        [TestCase(2, 2)]
        [TestCase(2, 3)]
        public void TestOfArray(int rows, int columns)
        {
            var sparseData = MatrixHelper.LoadSparse(rows, columns);

            var sparseA = sparseData.A;

            var array = new double[rows, columns];

            foreach (var a in sparseA.EnumerateIndexed())
            {
                array[a.Item1, a.Item2] = a.Item3;
            }

            var sparseB = SparseMatrix.OfArray(array);

            Assert.IsTrue(sparseA.Equals(sparseB));
        }

        [TestCase(2, 2)]
        [TestCase(2, 3)]
        public void TestOfJaggedArray(int rows, int columns)
        {
            var sparseData = MatrixHelper.LoadSparse(rows, columns);

            var sparseA = sparseData.A;

            var array = new double[rows][];

            for (int i = 0; i < rows; i++)
            {
                var values = new double[columns];

                for (int j = 0; j < columns; j++)
                {
                    values[j] = sparseA.At(i, j);
                }

                array[i] = values;
            }

            var sparseB = SparseMatrix.OfJaggedArray(array);

            Assert.IsTrue(sparseA.Equals(sparseB));
        }

        [Test]
        public void TestOfDiagonalArray()
        {
            int order = 3;

            var a = 1.0;
            var diag = Vector.Create(order, a);

            var A = SparseMatrix.OfDiagonalArray(diag);

            for (int i = 0; i < order; i++)
            {
                Assert.AreEqual(A.At(i, i), a);
            }
        }

        [Test]
        public void TestCreateIdentity()
        {
            int order = 3;

            var A = SparseMatrix.CreateIdentity(order);

            for (int i = 0; i < order; i++)
            {
                Assert.AreEqual(A.At(i, i), 1.0);
            }
        }

        #endregion
    }
}
