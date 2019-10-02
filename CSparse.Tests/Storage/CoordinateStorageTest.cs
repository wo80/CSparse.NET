using CSparse.Storage;
using NUnit.Framework;

namespace CSparse.Tests.Storage
{
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

            var C = new CoordinateStorage<double>(10, 10, -1);

            Assert.IsNull(C.RowIndices);
            Assert.IsNull(C.ColumnIndices);
            Assert.IsNull(C.Values);
        }
    }
}
