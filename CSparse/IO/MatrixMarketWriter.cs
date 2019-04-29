namespace CSparse.IO
{
    using System;
    using System.Globalization;
    using System.IO;
    using System.Numerics;
    using CSparse.Storage;

    /// <summary>
    /// NIST MatrixMarket Format Writer (http://math.nist.gov/MatrixMarket/)
    /// </summary>
    public static class MatrixMarketWriter
    {
        static readonly NumberFormatInfo Format = CultureInfo.InvariantCulture.NumberFormat;

        /// <summary>
        /// Write a matrix to file.
        /// </summary>
        /// <param name="filePath">The file path.</param>
        /// <param name="matrix">The matrix to write.</param>
        public static void WriteMatrix<T>(string filePath, Matrix<T> matrix)
            where T : struct, IEquatable<T>, IFormattable
        {
            using (var stream = File.Create(filePath))
            {
                using (var writer = new StreamWriter(stream))
                {
                    WriteMatrix(writer, matrix);
                }
            }
        }

        /// <summary>
        /// Write a matrix to file.
        /// </summary>
        /// <param name="stream">The stream to write to.</param>
        /// <param name="matrix">The matrix to write.</param>
        public static void WriteMatrix<T>(Stream stream, Matrix<T> matrix)
            where T : struct, IEquatable<T>, IFormattable
        {
            using (var writer = new StreamWriter(stream))
            {
                WriteMatrix(writer, matrix);
            }
        }

        /// <summary>
        /// Write a matrix to file.
        /// </summary>
        /// <param name="writer">The stream to write to.</param>
        /// <param name="matrix">The matrix to write.</param>
        public static void WriteMatrix<T>(StreamWriter writer, Matrix<T> matrix)
            where T : struct, IEquatable<T>, IFormattable
        {
            var complex = typeof(T) == typeof(Complex);
            var format = CreateValueFormatter<T>();

            var sparse = matrix as CompressedColumnStorage<T>;
            if (sparse != null)
            {
                writer.WriteLine("%%MatrixMarket matrix coordinate {0} general", complex ? "complex" : "real");
                writer.WriteLine("{0} {1} {2}", sparse.RowCount, sparse.ColumnCount, sparse.NonZerosCount);
                for (int column = 0; column < sparse.ColumnCount; column++)
                {
                    var endIndex = sparse.ColumnPointers[column + 1];
                    for (var j = sparse.ColumnPointers[column]; j < endIndex; j++)
                    {
                        writer.WriteLine("{0} {1} {2}", sparse.RowIndices[j] + 1, column + 1, format(sparse.Values[j]));
                    }
                }

                return;
            }
            
            var dense = matrix as DenseColumnMajorStorage<T>;
            if (dense != null)
            {
                writer.WriteLine("%%MatrixMarket matrix array {0} general", complex ? "complex" : "real");
                writer.WriteLine("{0} {1}", dense.RowCount, dense.ColumnCount);
                foreach (var value in dense.Values)
                {
                    writer.WriteLine(format(value));
                }

                return;
            }
        }
        
        static Func<T, string> CreateValueFormatter<T>() where T : IFormattable
        {
            if (typeof(T) == typeof(double))
            {
                return value => string.Format(Format, "{0:G14}", value);
            }

            if (typeof(T) == typeof(float))
            {
                return value => string.Format(Format, "{0:G7}", value);
            }

            if (typeof(T) == typeof(Complex))
            {
                return value =>
                {
                    var c = (Complex)(object)value;
                    return string.Format(Format, "{0:G14} {1:G14}", c.Real, c.Imaginary);
                };
            }

            throw new NotSupportedException();
        }
    }
}
