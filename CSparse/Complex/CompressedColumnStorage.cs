
namespace CSparse.Complex
{
    using System;

    [Obsolete("Use SparseMatrix class.")]
    public class CompressedColumnStorage : SparseMatrix
    {
        /// <summary>
        /// Initializes a new instance of the CompressedColumnStorage class.
        /// </summary>
        public CompressedColumnStorage(int rowCount, int columnCount)
            : base(rowCount, columnCount)
        {
        }

        /// <summary>
        /// Initializes a new instance of the CompressedColumnStorage class.
        /// </summary>
        public CompressedColumnStorage(int rowCount, int columnCount, int valueCount)
            : base(rowCount, columnCount, valueCount)
        {
        }
    }
}
