
namespace CSparse.Tests.IO
{
    using CSparse.IO;
    using NUnit.Framework;
    using System.IO;

    public class TestMatrixMarketReader
    {
        [Test]
        public void TestSymmetric()
        {
            var stream = ResourceLoader.GetStream("LFAT5.mtx", "Double");

            using (var reader = new StreamReader(stream))
            {
                var A = MatrixMarketReader.ReadStorage<double>(reader);

                Assert.AreEqual(A.RowCount, 14);
                Assert.AreEqual(A.ColumnCount, 14);

                // Symmetric mtx file has 30 entries -> auto expand = 2 * 30 - 14 = 46.
                Assert.AreEqual(A.NonZerosCount, 46);
            }
        }

        [Test]
        public void TestPatternSymmetric()
        {
            var stream = ResourceLoader.GetStream("bcspwr01.mtx", "Double");

            using (var reader = new StreamReader(stream))
            {
                var A = MatrixMarketReader.ReadStorage<double>(reader, false);

                Assert.AreEqual(A.RowCount, 39);
                Assert.AreEqual(A.ColumnCount, 39);

                // Symmetric mtx file has 85 entries, no auto expand.
                Assert.AreEqual(A.NonZerosCount, 85);
            }
        }
    }
}
