
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

            Assert.IsNotNull(A.RowIndices);
            Assert.IsNotNull(A.ColumnIndices);
            Assert.IsNotNull(A.Values);

            A.At(1, 2, 3.0);

            var B = new CoordinateStorage<double>(10, 10, 0);

            Assert.IsNotNull(B.RowIndices);
            Assert.IsNotNull(B.ColumnIndices);
            Assert.IsNotNull(B.Values);

            B.At(1, 2, 3.0);

            Assert.Throws<ArgumentOutOfRangeException>(() => new CoordinateStorage<double>(-1,  1,  1));
            Assert.Throws<ArgumentOutOfRangeException>(() => new CoordinateStorage<double>( 1, -1,  1));
            Assert.Throws<ArgumentOutOfRangeException>(() => new CoordinateStorage<double>( 1,  1, -1));
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

                Assert.IsTrue(column < 1);
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

            Assert.IsTrue(AT.Equals(sparseB));

            cooT = coo.Transpose(false);

            var sparseC = SparseMatrix.OfIndexed(cooT);

            Assert.IsTrue(AT.Equals(sparseC));
            Assert.IsTrue(ReferenceEquals(coo.RowIndices, cooT.ColumnIndices));
            Assert.IsTrue(ReferenceEquals(coo.ColumnIndices, cooT.RowIndices));
            Assert.IsTrue(ReferenceEquals(coo.Values, cooT.Values));

        }
    }
}
