
namespace CSparse.Complex
{
    using System;

    [Obsolete("Use SparseMatrix class.")]
    public class CompressedColumnStorage : SparseMatrix
    {
        /// <inheritdoc />
        public CompressedColumnStorage(int rowCount, int columnCount)
            : base(rowCount, columnCount)
        {
        }

        /// <inheritdoc />
        public CompressedColumnStorage(int rowCount, int columnCount, int valueCount)
            : base(rowCount, columnCount, valueCount)
        {
        }
    }
}
