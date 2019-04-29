namespace CSparse.Tests
{
    using CSparse.Storage;
    using System;

    class SparseTestData<T>
        where T : struct, IEquatable<T>, IFormattable
    {
        public int RowCount;

        public int ColumnCount;

        public CompressedColumnStorage<T> A;

        public CompressedColumnStorage<T> B;

        public T[] x;

        public T[] y;

        public CompressedColumnStorage<T> AT;

        public CompressedColumnStorage<T> BT;

        public CompressedColumnStorage<T> ApB;

        public CompressedColumnStorage<T> AmBT;

        public CompressedColumnStorage<T> ATmB;

        public T[] Ax;

        public T[] ATy;

        public T[] xTBT;
    }
}
