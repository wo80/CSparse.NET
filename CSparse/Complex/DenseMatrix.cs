
namespace CSparse.Complex
{
    using CSparse.Properties;
    using CSparse.Storage;
    using System;
    using System.Diagnostics;
    using System.Numerics;

    /// <summary>
    /// Dense matrix stored in column major order.
    /// </summary>
    [DebuggerDisplay("DenseMatrix {RowCount}x{ColumnCount}")]
    [Serializable]
    public class DenseMatrix : DenseColumnMajorStorage<Complex>
    {
        /// <summary>
        /// Initializes a new instance of the DenseMatrix class.
        /// </summary>
        public DenseMatrix(int rows, int columns)
            : this(rows, columns, new Complex[rows * columns])
        {
        }

        /// <summary>
        /// Initializes a new instance of the DenseMatrix class.
        /// </summary>
        public DenseMatrix(int rows, int columns, Complex[] values)
            : base(rows, columns, values)
        {
        }

        /// <inheritdoc />
        public override double L1Norm()
        {
            double sum, norm = 0.0;

            for (var j = 0; j < columnCount; j++)
            {
                sum = 0.0;
                for (var i = 0; i < rowCount; i++)
                {
                    sum += values[(j * rowCount) + i].Magnitude;
                }
                norm = Math.Max(norm, sum);
            }

            return norm;
        }
        
        /// <inheritdoc />
        public override double InfinityNorm()
        {
            var r = new double[rowCount];

            for (var j = 0; j < columnCount; j++)
            {
                for (var i = 0; i < rowCount; i++)
                {
                    r[i] += values[(j * rowCount) + i].Magnitude;
                }
            }

            double norm = r[0];

            for (int i = 1; i < rowCount; i++)
            {
                if (r[i] > norm)
                {
                    norm = r[i];
                }
            }

            return norm;
        }

        /// <inheritdoc />
        public override double FrobeniusNorm()
        {
            double sum = 0.0, norm = 0.0;

            int length = rowCount * columnCount;

            for (int i = 0; i < length; i++)
            {
                sum = values[i].Magnitude;
                norm += sum * sum;
            }

            return Math.Sqrt(norm);
        }

        /// <inheritdoc />
        public override object Clone()
        {
            var values = (Complex[])this.values.Clone();

            return new DenseMatrix(rowCount, columnCount, values);
        }

        /// <inheritdoc />
        public override void Multiply(Complex[] x, Complex[] y)
        {
            var A = values;

            int rows = rowCount;
            int cols = columnCount;

            for (int i = 0; i < rows; i++)
            {
                Complex sum = Complex.Zero;

                for (int j = 0; j < cols; j++)
                {
                    sum += A[(j * rows) + i] * x[j];
                }

                y[i] = sum;
            }
        }

        /// <inheritdoc />
        public override void Multiply(Complex alpha, Complex[] x, Complex beta, Complex[] y)
        {
            var A = values;

            int rows = rowCount;
            int cols = columnCount;

            for (int i = 0; i < rows; i++)
            {
                Complex sum = Complex.Zero;

                for (int j = 0; j < cols; j++)
                {
                    sum += A[(j * rows) + i] * x[j];
                }

                y[i] = beta * y[i] + alpha * sum;
            }
        }

        /// <inheritdoc />
        public override void TransposeMultiply(Complex[] x, Complex[] y)
        {
            var A = values;

            int rows = rowCount;
            int cols = columnCount;

            for (int j = 0; j < cols; j++)
            {
                int row = j * rows;

                Complex sum = Complex.Zero;

                for (int i = 0; i < rows; i++)
                {
                    sum += A[row + i] * x[i];
                }

                y[j] = sum;
            }
        }

        /// <inheritdoc />
        public override void TransposeMultiply(Complex alpha, Complex[] x, Complex beta, Complex[] y)
        {
            var A = values;

            int rows = rowCount;
            int cols = columnCount;

            for (int j = 0; j < cols; j++)
            {
                int row = j * rows;

                Complex sum = Complex.Zero;

                for (int i = 0; i < rows; i++)
                {
                    sum += A[row + i] * x[i];
                }

                y[j] = beta * y[j] + alpha * sum;
            }
        }

        /// <inheritdoc />
        public override void Add(DenseColumnMajorStorage<Complex> other, DenseColumnMajorStorage<Complex> result)
        {
            int m = rowCount;
            int n = columnCount;

            // check inputs
            if (m != other.RowCount || n != other.ColumnCount)
            {
                throw new ArgumentException();
            }

            var target = result.Values;

            var a = this.values;
            var b = other.Values;

            int length = m * n;

            for (int i = 0; i < length; i++)
            {
                target[i] = a[i] + b[i];
            }
        }

        /// <inheritdoc />
        public override void Multiply(DenseColumnMajorStorage<Complex> other, DenseColumnMajorStorage<Complex> result)
        {
            var A = values;
            var B = other.Values;
            var C = result.Values;

            int m = rowCount; // rows of matrix A
            int n = other.ColumnCount;
            int o = columnCount;

            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    Complex sum = Complex.Zero;

                    for (int k = 0; k < o; ++k)
                    {
                        sum += A[(k * m) + i] * B[(j * o) + k];
                    }

                    C[(j * m) + i] += sum;
                }
            }
        }

        /// <inheritdoc />
        public override void PointwiseMultiply(DenseColumnMajorStorage<Complex> other, DenseColumnMajorStorage<Complex> result)
        {
            if (RowCount != other.RowCount || ColumnCount != other.ColumnCount)
            {
                throw new ArgumentException(Resources.MatrixDimensions);
            }

            if (RowCount != result.RowCount || ColumnCount != result.ColumnCount)
            {
                throw new ArgumentException(Resources.MatrixDimensions);
            }

            var x = values;
            var y = other.Values;

            var target = result.Values;

            int length = target.Length;

            for (int i = 0; i < length; i++)
            {
                target[i] = x[i] * y[i];
            }
        }

        /// <inheritdoc />
        public override bool Equals(Matrix<Complex> other, double tolerance)
        {
            if (rowCount != other.RowCount || columnCount != other.ColumnCount)
            {
                return false;
            }

            var dense = other as DenseColumnMajorStorage<Complex>;

            if (dense == null)
            {
                return false;
            }

            int length = rowCount * columnCount;

            var otherValues = dense.Values;

            for (int i = 0; i < length; i++)
            {
                if (Complex.Abs(values[i] - otherValues[i]) > tolerance)
                {
                    return false;
                }
            }

            return true;
        }
    }
}
