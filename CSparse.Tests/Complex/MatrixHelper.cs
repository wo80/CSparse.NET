namespace CSparse.Tests.Complex
{
    using CSparse.Complex;
    using CSparse.Storage;
    using System.Collections.Generic;
    using System.IO;
    using System.Linq;
    using System.Numerics;

    class MatrixHelper
    {
        private static Dictionary<string, DenseTestData<Complex>> dense = new Dictionary<string, DenseTestData<Complex>>();

        private static Dictionary<string, SparseTestData<Complex>> sparse = new Dictionary<string, SparseTestData<Complex>>();

        public static SparseTestData<Complex> LoadSparse(int rows, int columns)
        {
            string resource = string.Format("test-data-dense-{0}x{1}.txt", rows, columns);

            SparseTestData<Complex> data;

            if (!sparse.TryGetValue(resource, out data))
            {
                var dense = LoadDense(rows, columns);

                data = DenseToSparse(dense);

                sparse.Add(resource, data);
            }

            return data;
        }

        public static DenseTestData<Complex> LoadDense(int rows, int columns)
        {
            string resource = string.Format("test-data-dense-{0}x{1}.txt", rows, columns);

            DenseTestData<Complex> data;

            if (!dense.TryGetValue(resource, out data))
            {
                var stream = ResourceLoader.GetStream(resource, "Double");

                data = ReadDenseTestData(stream);

                dense.Add(resource, data);
            }

            return data;
        }

        private static DenseTestData<Complex> ReadDenseTestData(Stream stream)
        {
            var data = Tests.Double.DenseTestDataReader.Read(stream);

            var result = new DenseTestData<Complex>();

            result.A = ToComplex(data.A);
            result.B = ToComplex(data.B);
            result.x = ToComplex(data.x);
            result.y = ToComplex(data.y);
            result.AT = ToComplex(data.AT);
            result.BT = ToComplex(data.BT);
            result.ApB = ToComplex(data.ApB);
            result.AmBT = ToComplex(data.AmBT);
            result.ATmB = ToComplex(data.ATmB);
            result.Ax = ToComplex(data.Ax);
            result.ATy = ToComplex(data.ATy);
            result.xTBT = ToComplex(data.xTBT);

            return result;
        }

        private static Complex[] ToComplex(double[] vec)
        {
            int length = vec.Length;

            var result = new Complex[length];

            for (int i = 0; i < length; i++)
            {
                result[i] = vec[i];
            }

            return result;
        }

        private static DenseColumnMajorStorage<Complex> ToComplex(DenseColumnMajorStorage<double> matrix)
        {
            var result = new DenseMatrix(matrix.RowCount, matrix.ColumnCount);

            int length = matrix.RowCount * matrix.ColumnCount;

            for (int i = 0; i < length; i++)
            {
                result.Values[i] = matrix.Values[i];
            }

            return result;
        }

        private static SparseTestData<Complex> DenseToSparse(DenseTestData<Complex> dense)
        {
            var data = new SparseTestData<Complex>()
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

        private static CompressedColumnStorage<Complex> DenseToSparse(DenseColumnMajorStorage<Complex> dense)
        {
            var cs = Converter.FromColumnMajorArray(dense.Values, dense.RowCount, dense.ColumnCount);

            return Converter.ToCompressedColumnStorage(cs);
        }
    }
}
