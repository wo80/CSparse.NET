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
        /// Initializes a new instance of the <see cref="DenseMatrix"/> class.
        /// </summary>
        /// <param name="size">The size of the square matrix.</param>
        public DenseMatrix(int size)
            : base(size, size)
        {
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="DenseMatrix"/> class.
        /// </summary>
        public DenseMatrix(int rows, int columns)
            : base(rows, columns)
        {
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="DenseMatrix"/> class.
        /// </summary>
        public DenseMatrix(int rows, int columns, Complex[] values)
            : base(rows, columns, values)
        {
        }

        /// <inheritdoc />
        public override double L1Norm()
        {
            double sum, norm = 0.0;

            for (var j = 0; j < columns; j++)
            {
                sum = 0.0;
                for (var i = 0; i < rows; i++)
                {
                    sum += Values[(j * rows) + i].Magnitude;
                }
                norm = Math.Max(norm, sum);
            }

            return norm;
        }
        
        /// <inheritdoc />
        public override double InfinityNorm()
        {
            var r = new double[rows];

            for (var j = 0; j < columns; j++)
            {
                for (var i = 0; i < rows; i++)
                {
                    r[i] += Values[(j * rows) + i].Magnitude;
                }
            }

            double norm = r[0];

            for (int i = 1; i < rows; i++)
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

            int length = rows * columns;

            for (int i = 0; i < length; i++)
            {
                sum = Values[i].Magnitude;
                norm += sum * sum;
            }

            return Math.Sqrt(norm);
        }

        /// <inheritdoc />
        public override DenseColumnMajorStorage<Complex> Clone()
        {
            var values = (Complex[])this.Values.Clone();

            return new DenseMatrix(rows, columns, values);
        }

        /// <inheritdoc />
        public override void Multiply(ReadOnlySpan<Complex> x, Span<Complex> y)
        {
            var A = Values;

            int rows = base.rows;
            int cols = columns;

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
        public override void Multiply(Complex alpha, ReadOnlySpan<Complex> x, Complex beta, Span<Complex> y)
        {
            var A = Values;

            int rows = base.rows;
            int cols = columns;

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
        public override void TransposeMultiply(ReadOnlySpan<Complex> x, Span<Complex> y)
        {
            var A = Values;

            int rows = base.rows;
            int cols = columns;

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
        public override void TransposeMultiply(Complex alpha, ReadOnlySpan<Complex> x, Complex beta, Span<Complex> y)
        {
            var A = Values;

            int rows = base.rows;
            int cols = columns;

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
            int m = rows;
            int n = columns;

            // check inputs
            if (m != other.RowCount || n != other.ColumnCount)
            {
                throw new ArgumentException();
            }

            var target = result.Values;

            var a = this.Values;
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
            var A = Values;
            var B = other.Values;
            var C = result.Values;

            int m = rows; // rows of matrix A
            int n = other.ColumnCount;
            int o = columns;

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
                throw new ArgumentException(Resources.MatrixDimensions, nameof(other));
            }

            if (RowCount != result.RowCount || ColumnCount != result.ColumnCount)
            {
                throw new ArgumentException(Resources.MatrixDimensions, nameof(result));
            }

            var x = Values;
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
            if (rows != other.RowCount || columns != other.ColumnCount)
            {
                return false;
            }

            var dense = other as DenseColumnMajorStorage<Complex>;

            if (dense == null)
            {
                return false;
            }

            int length = rows * columns;

            var otherValues = dense.Values;

            for (int i = 0; i < length; i++)
            {
                if (Complex.Abs(Values[i] - otherValues[i]) > tolerance)
                {
                    return false;
                }
            }

            return true;
        }
    }
}
