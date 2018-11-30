
namespace CSparse.Tests.Double
{
    using CSparse.Double;
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

        #endregion

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
        public void TestOfColumnMajor(int rows, int columns)
        {
            var denseData = MatrixHelper.LoadDense(rows, columns);
            var denseA = denseData.A;

            var sparseData = MatrixHelper.LoadSparse(rows, columns);
            var sparseA = sparseData.A;

            var sparseB = SparseMatrix.OfColumnMajor(rows, columns, denseA.Values);

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
