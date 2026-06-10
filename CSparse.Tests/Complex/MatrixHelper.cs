namespace CSparse.Tests.Complex
{
    using CSparse.Complex;
    using CSparse.Storage;
    using System.Collections.Generic;
    using System.IO;
    using System.Numerics;

    class MatrixHelper
    {
        private static readonly Dictionary<string, DenseTestData<Complex>> dense = [];
        private static readonly Dictionary<string, SparseTestData<Complex>> sparse = [];

        public static SparseTestData<Complex> LoadSparse(int rows, int columns)
        {
            string resource = string.Format("test-data-dense-{0}x{1}.txt", rows, columns);

            if (!sparse.TryGetValue(resource, out SparseTestData<Complex> data))
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

            if (!dense.TryGetValue(resource, out DenseTestData<Complex> data))
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

            return new DenseTestData<Complex>
            {
                A = ToComplex(data.A),
                B = ToComplex(data.B),
                x = ToComplex(data.x),
                y = ToComplex(data.y),
                AT = ToComplex(data.AT),
                BT = ToComplex(data.BT),
                ApB = ToComplex(data.ApB),
                AmBT = ToComplex(data.AmBT),
                ATmB = ToComplex(data.ATmB),
                Ax = ToComplex(data.Ax),
                ATy = ToComplex(data.ATy),
                xTBT = ToComplex(data.xTBT)
            };
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
            return new SparseTestData<Complex>
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

        private static CompressedColumnStorage<Complex> DenseToSparse(DenseColumnMajorStorage<Complex> dense)
        {
            return CompressedColumnStorage<Complex>.OfColumnMajor(dense.RowCount, dense.ColumnCount, dense.Values);
        }
    }
}
