
namespace CSparse.Tests
{
    using CSparse.Storage;
    using System;

    class DenseTestData<T>
        where T : struct, IEquatable<T>, IFormattable
    {
        public int RowCount;

        public int ColumnCount;

        public DenseColumnMajorStorage<T> A;

        public DenseColumnMajorStorage<T> B;

        public T[] x;

        public T[] y;

        public DenseColumnMajorStorage<T> AT;

        public DenseColumnMajorStorage<T> BT;

        public DenseColumnMajorStorage<T> ApB;

        public DenseColumnMajorStorage<T> AmBT;

        public DenseColumnMajorStorage<T> ATmB;

        public T[] Ax;

        public T[] ATy;

        public T[] xTBT;
    }
}
