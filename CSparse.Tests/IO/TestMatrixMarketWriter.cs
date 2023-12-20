
namespace CSparse.Tests.IO
{
    using System.IO;
    using CSparse.Double;
    using CSparse.IO;
    using NUnit.Framework;

    public class TestMatrixMarketWriter
    {
        private static readonly double[] DenseA = new[]
        {
            1.0, 0.0, 0.5,
            0.0, 2.0, 0.1,
            0.5, 0.1, 3.0,
        };

        [Test]
        public void TestWriteSparse()
        {
            var A = SparseMatrix.OfColumnMajor(3, 3, DenseA);

            using var stream = new MemoryStream(256);
            using var writer = new StreamWriter(stream, leaveOpen: true);

            MatrixMarketWriter.WriteMatrix(writer, A);

            stream.Position = 0;

            var B = MatrixMarketReader.ReadMatrix<double>(stream);

            Assert.That(A.Equals(B), Is.True);
        }

        [Test]
        public void TestWriteSparseSymmetric()
        {
            var A = SparseMatrix.OfColumnMajor(3, 3, DenseA);

            using var stream = new MemoryStream(256);
            using var writer = new StreamWriter(stream, leaveOpen: true);

            MatrixMarketWriter.WriteMatrix(writer, A, true);

            stream.Position = 0;

            var B = MatrixMarketReader.ReadMatrix<double>(stream);

            Assert.That(A.Equals(B), Is.True);
        }

        [Test]
        public void TestWriteDense()
        {
            var A = DenseMatrix.OfColumnMajor(3, 3, DenseA);

            using var stream = new MemoryStream(256);
            using var writer = new StreamWriter(stream, leaveOpen: true);

            MatrixMarketWriter.WriteMatrix(writer, A);

            // MatrixMarketReader doesn't support dense (array) storage,
            // so we just check that some bytes were written.
            Assert.That(stream.Position, Is.GreaterThan(0));
        }
    }
}
