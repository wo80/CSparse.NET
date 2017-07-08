
namespace CSparse.Tests.Complex
{
    using CSparse.Complex;
    using System;

    class MatrixHelper
    {
        public static SparseMatrix Load(int rows, int columns)
        {
            if (rows == 0 || columns == 0)
            {
                return new SparseMatrix(rows, columns, 0);
            }

            throw new NotImplementedException();
        }
    }
}
