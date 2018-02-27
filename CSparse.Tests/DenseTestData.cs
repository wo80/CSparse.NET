
namespace CSparse.Tests
{
    using System;

    class DenseTestData<T>
        where T : struct, IEquatable<T>, IFormattable
    {
        public int RowCount;

        public int ColumnCount;

        public T[,] A;

        public T[,] B;

        public T[] x;

        public T[] y;

        public T[,] AT;

        public T[,] BT;

        public T[,] ApB;

        public T[,] AmBT;

        public T[,] ATmB;

        public T[] Ax;

        public T[] ATy;

        public T[] xTBT;
    }
}
