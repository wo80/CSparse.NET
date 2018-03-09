
namespace CSparse.Tests.Double
{
    using CSparse.Storage;
    using System.Collections.Generic;

    class MatrixHelper
    {
        private static Dictionary<string, DenseTestData<double>> dense = new Dictionary<string, DenseTestData<double>>();

        private static Dictionary<string, SparseTestData<double>> sparse = new Dictionary<string, SparseTestData<double>>();

        public static SparseTestData<double> LoadSparse(int rows, int columns)
        {
            string resource = string.Format("test-data-dense-{0}x{1}.txt", rows, columns);

            SparseTestData<double> data;

            if (!sparse.TryGetValue(resource, out data))
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

            DenseTestData<double> data;

            if (!dense.TryGetValue(resource, out data))
            {
                var stream = ResourceLoader.GetStream(resource, "Double");

                data = DenseTestDataReader.Read(stream);

                dense.Add(resource, data);
            }

            return data;
        }

        private static SparseTestData<double> DenseToSparse(DenseTestData<double> dense)
        {
            var data = new SparseTestData<double>()
            {
                RowCount = dense.RowCount,
                ColumnCount = dense.ColumnCount
            };

            data.A = DenseToSparse(dense.A);
            data.B = DenseToSparse(dense.B);
            data.x = dense.x;
            data.y = dense.y;
            data.AT = DenseToSparse(dense.AT);
            data.BT = DenseToSparse(dense.BT);
            data.ApB = DenseToSparse(dense.ApB);
            data.AmBT = DenseToSparse(dense.AmBT);
            data.ATmB = DenseToSparse(dense.ATmB);
            data.Ax = dense.Ax;
            data.ATy = dense.ATy;
            data.xTBT = dense.xTBT;

            return data;
        }

        private static CompressedColumnStorage<double> DenseToSparse(DenseColumnMajorStorage<double> dense)
        {
            var cs = Converter.FromColumnMajorArray(dense.Values, dense.RowCount, dense.ColumnCount);

            return Converter.ToCompressedColumnStorage(cs);
        }
    }
}
