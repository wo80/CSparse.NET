
namespace CSparse.Tests.Storage
{
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
    }
}
