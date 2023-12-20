
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

                Assert.That(14, Is.EqualTo(A.RowCount));
                Assert.That(14, Is.EqualTo(A.ColumnCount));

                // Symmetric mtx file has 30 entries -> auto expand = 2 * 30 - 14 = 46.
                Assert.That(46, Is.EqualTo(A.NonZerosCount));
            }
        }

        [Test]
        public void TestPatternSymmetric()
        {
            var stream = ResourceLoader.GetStream("bcspwr01.mtx", "Double");

            using (var reader = new StreamReader(stream))
            {
                var A = MatrixMarketReader.ReadStorage<double>(reader, false);

                Assert.That(39, Is.EqualTo(A.RowCount));
                Assert.That(39, Is.EqualTo(A.ColumnCount));

                // Symmetric mtx file has 85 entries, no auto expand.
                Assert.That(85, Is.EqualTo(A.NonZerosCount));
            }
        }
    }
}
