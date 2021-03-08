namespace CSparse.Double
{
    using CSparse.Properties;
    using CSparse.Storage;
    using System;
    using System.Diagnostics;
    using System.Threading.Tasks;

    /// <summary>
    /// Dense matrix stored in column major order.
    /// </summary>
    [DebuggerDisplay("DenseMatrix {RowCount}x{ColumnCount}")]
    [Serializable]
    public class DenseMatrix : DenseColumnMajorStorage<Double>
    {
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
        public DenseMatrix(int rows, int columns, double[] values)
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
                    sum += Math.Abs(Values[(j * rows) + i]);
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
                    r[i] += Math.Abs(Values[(j * rows) + i]);
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
                sum = Math.Abs(Values[i]);
                norm += sum * sum;
            }

            return Math.Sqrt(norm);
        }

        /// <inheritdoc />
        public override DenseColumnMajorStorage<double> Clone()
        {
            var values = (double[])this.Values.Clone();

            return new DenseMatrix(rows, columns, values);
        }

        /// <inheritdoc />
        public override void Multiply(double[] x, double[] y)
        {
            var A = Values;

            int rows = base.rows;
            int cols = columns;

            for (int i = 0; i < rows; i++)
            {
                double sum = 0.0;

                for (int j = 0; j < cols; j++)
                {
                    sum += A[(j * rows) + i] * x[j];
                }

                y[i] = sum;
            }
        }

        /// <inheritdoc />
        public override void Multiply(double alpha, double[] x, double beta, double[] y)
        {
            var A = Values;

            int rows = base.rows;
            int cols = columns;

            for (int i = 0; i < rows; i++)
            {
                double sum = 0.0;

                for (int j = 0; j < cols; j++)
                {
                    sum += A[(j * rows) + i] * x[j];
                }

                y[i] = beta * y[i] + alpha * sum;
            }
        }

        /// <inheritdoc />
        public override void TransposeMultiply(double[] x, double[] y)
        {
            var A = Values;

            int rows = base.rows;
            int cols = columns;
            
            for (int j = 0; j < cols; j++)
            {
                int col = j * rows;

                double sum = 0.0;

                for (int i = 0; i < rows; i++)
                {
                    sum += A[col + i] * x[i];
                }

                y[j] = sum;
            }
        }

        /// <inheritdoc />
        public override void TransposeMultiply(double alpha, double[] x, double beta, double[] y)
        {
            var A = Values;

            int rows = base.rows;
            int cols = columns;

            for (int j = 0; j < cols; j++)
            {
                y[j] = beta * y[j];
            }
            
            for (int j = 0; j < cols; j++)
            {
                int col = j * rows;

                double sum = 0.0;

                for (int i = 0; i < rows; i++)
                {
                    sum += A[col + i] * x[i];
                }

                y[j] = beta * y[j] +  alpha * sum;
            }
        }

        /// <inheritdoc />
        public override void Add(DenseColumnMajorStorage<double> other, DenseColumnMajorStorage<double> result)
        {
            int rows = this.rows;
            int columns = this.columns;

            // check inputs
            if (rows != other.RowCount || columns != other.ColumnCount)
            {
                throw new ArgumentException();
            }

            var target = result.Values;

            var a = this.Values;
            var b = other.Values;

            int length = rows * columns;

            for (int i = 0; i < length; i++)
            {
                target[i] = a[i] + b[i];
            }
        }

        /// <inheritdoc />
        public override void Multiply(DenseColumnMajorStorage<double> other, DenseColumnMajorStorage<double> result)
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
                    double sum = 0.0;

                    for (int k = 0; k < o; ++k)
                    {
                        sum += A[(k * m) + i] * B[(j * o) + k];
                    }

                    C[(j * m) + i] += sum;
                }
            }
        }

        /// <inheritdoc />
        public override void ParallelMultiply(DenseColumnMajorStorage<double> other, DenseColumnMajorStorage<double> result, ParallelOptions options = null)
        {
            var A = Values;
            var B = other.Values;
            var C = result.Values;

            int m = rows; // rows of matrix A
            int n = other.ColumnCount;
            int o = columns;

            int processorCount = Environment.ProcessorCount;

            // Allow for at least 2 threads with 4 rows each
            if (m < 2 * 4 || n <= 0 || processorCount < 2 || (long)m * n * o < 1000000L)
            {
                Multiply(other, result);
                return;
            }

            if (options == null)
            {
                options = new ParallelOptions() { MaxDegreeOfParallelism = processorCount };
            }
            else if (options.MaxDegreeOfParallelism < 0 || options.MaxDegreeOfParallelism > processorCount)
            {
                options.MaxDegreeOfParallelism = processorCount;
            }

            var nblocks = Math.Min(options.MaxDegreeOfParallelism, m);
            var starts = new int[nblocks + 1];
            for (var i = 0; i <= nblocks; i++)
            {
                starts[i] = i * m / nblocks;
            }

            Parallel.For(0, nblocks, options, index =>
            {
                for (int i = starts[index]; i < starts[index + 1]; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        double sum = 0.0;

                        for (int k = 0; k < o; ++k)
                        {
                            sum += A[(k * m) + i] * B[(j * o) + k];
                        }

                        C[(j * m) + i] += sum;
                    }
                }
            });
        }

        /// <inheritdoc />
        public override void PointwiseMultiply(DenseColumnMajorStorage<double> other, DenseColumnMajorStorage<double> result)
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
        public override bool Equals(Matrix<double> other, double tolerance)
        {
            if (rows != other.RowCount || columns != other.ColumnCount)
            {
                return false;
            }

            var dense = other as DenseColumnMajorStorage<double>;

            if (dense == null)
            {
                return false;
            }

            int length = rows * columns;

            var otherValues = dense.Values;

            for (int i = 0; i < length; i++)
            {
                if (Math.Abs(Values[i] - otherValues[i]) > tolerance)
                {
                    return false;
                }
            }

            return true;
        }
    }
}
