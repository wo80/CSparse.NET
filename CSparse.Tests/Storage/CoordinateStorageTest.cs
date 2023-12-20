
namespace CSparse.Tests.Storage
{
    using CSparse.Double;
    using CSparse.Storage;
    using NUnit.Framework;
    using System;

    public class CoordinateStorageTest
    {
        [Test]
        public void TestConstructor()
        {
            var A = new CoordinateStorage<double>(10, 10, 1);

            Assert.That(A.RowIndices, Is.Not.Null);
            Assert.That(A.ColumnIndices, Is.Not.Null);
            Assert.That(A.Values, Is.Not.Null);

            A.At(1, 2, 3.0);

            var B = new CoordinateStorage<double>(10, 10, 0);

            Assert.That(B.RowIndices, Is.Not.Null);
            Assert.That(B.ColumnIndices, Is.Not.Null);
            Assert.That(B.Values, Is.Not.Null);

            B.At(1, 2, 3.0);

            Assert.Throws<ArgumentOutOfRangeException>(() => new CoordinateStorage<double>(-1,  1,  1));
            Assert.Throws<ArgumentOutOfRangeException>(() => new CoordinateStorage<double>( 1, -1,  1));
            Assert.Throws<ArgumentOutOfRangeException>(() => new CoordinateStorage<double>( 1,  1, -1));

            // With storage arrays.

            Assert.Throws<ArgumentException>(() => new CoordinateStorage<double>(2, 2, 2, new int[2], new int[2], new double[1]));
        }

        [TestCase(2, 2)]
        [TestCase(2, 3)]
        public void TestKeep(int rows, int columns)
        {
            var sparseData = Tests.Double.MatrixHelper.LoadSparse(rows, columns);

            var sparseA = sparseData.A;

            var coo = new CoordinateStorage<double>(rows, columns, sparseA.NonZerosCount);

            foreach (var a in sparseA.EnumerateIndexed())
            {
                coo.At(a.Item1, a.Item2, a.Item3);
            }

            // Only keep the first column.
            coo.Keep((i, j, a) => j < 1);

            var sparseB = SparseMatrix.OfIndexed(coo);

            foreach (var a in sparseB.EnumerateIndexed())
            {
                int column = a.Item2;

                Assert.That(column < 1, Is.True);
            }
        }

        [TestCase(2, 2)]
        [TestCase(2, 3)]
        public void TestTranspose(int rows, int columns)
        {
            var sparseData = Tests.Double.MatrixHelper.LoadSparse(rows, columns);

            var sparseA = sparseData.A;
            var AT = sparseA.Transpose();

            var coo = new CoordinateStorage<double>(rows, columns, sparseA.NonZerosCount);

            foreach (var a in sparseA.EnumerateIndexed())
            {
                coo.At(a.Item1, a.Item2, a.Item3);
            }

            var cooT = coo.Transpose(true);

            var sparseB = SparseMatrix.OfIndexed(cooT);

            Assert.That(AT.Equals(sparseB), Is.True);

            cooT = coo.Transpose(false);

            var sparseC = SparseMatrix.OfIndexed(cooT);

            Assert.That(AT.Equals(sparseC), Is.True);
            Assert.That(ReferenceEquals(coo.RowIndices, cooT.ColumnIndices), Is.True);
            Assert.That(ReferenceEquals(coo.ColumnIndices, cooT.RowIndices), Is.True);
            Assert.That(ReferenceEquals(coo.Values, cooT.Values), Is.True);

        }
    }
}
