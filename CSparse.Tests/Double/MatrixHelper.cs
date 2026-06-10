namespace CSparse.Tests.Double
{
    using CSparse.Storage;
    using System.Collections.Generic;

    class MatrixHelper
    {
        private static readonly Dictionary<string, DenseTestData<double>> dense = [];
        private static readonly Dictionary<string, SparseTestData<double>> sparse = [];

        public static SparseTestData<double> LoadSparse(int rows, int columns)
        {
            string resource = string.Format("test-data-dense-{0}x{1}.txt", rows, columns);

            if (!sparse.TryGetValue(resource, out SparseTestData<double> data))
            {
                var dense = LoadDense(rows, columns);

                data = DenseToSparse(dense);

                sparse.Add(resource, data);
            }

            return data;
        }

        public static DenseTestData<double> LoadDense(int rows, int columns)
        {
            string resource = string.Format("test-data-dense-{0}x{1}.txt", rows, columns);

            if (!dense.TryGetValue(resource, out DenseTestData<double> data))
            {
                var stream = ResourceLoader.GetStream(resource, "Double");

                data = DenseTestDataReader.Read(stream);

                dense.Add(resource, data);
            }

            return data;
        }

        private static SparseTestData<double> DenseToSparse(DenseTestData<double> dense)
        {
            return new SparseTestData<double>
            {
                RowCount = dense.RowCount,
                ColumnCount = dense.ColumnCount,
                A = DenseToSparse(dense.A),
                B = DenseToSparse(dense.B),
                x = dense.x,
                y = dense.y,
                AT = DenseToSparse(dense.AT),
                BT = DenseToSparse(dense.BT),
                ApB = DenseToSparse(dense.ApB),
                AmBT = DenseToSparse(dense.AmBT),
                ATmB = DenseToSparse(dense.ATmB),
                Ax = dense.Ax,
                ATy = dense.ATy,
                xTBT = dense.xTBT
            };
        }

        private static CompressedColumnStorage<double> DenseToSparse(DenseColumnMajorStorage<double> dense)
        {
            return CompressedColumnStorage<double>.OfColumnMajor(dense.RowCount, dense.ColumnCount, dense.Values);
        }
    }
}
