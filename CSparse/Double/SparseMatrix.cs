namespace CSparse.Double
{
    using CSparse.Properties;
    using CSparse.Storage;
    using System;
    using System.Diagnostics;
    using System.Threading.Tasks;

    /// <inheritdoc />
    [DebuggerDisplay("SparseMatrix {RowCount}x{ColumnCount}-Double {NonZerosCount}-NonZero")]
    [Serializable]
    public class SparseMatrix : CompressedColumnStorage<double>
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="SparseMatrix"/> class.
        /// </summary>
        public SparseMatrix(int rowCount, int columnCount)
            : base(rowCount, columnCount)
        {
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="SparseMatrix"/> class.
        /// </summary>
        public SparseMatrix(int rowCount, int columnCount, int valueCount)
            : base(rowCount, columnCount, valueCount)
        {
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="SparseMatrix"/> class.
        /// </summary>
        public SparseMatrix(int rowCount, int columnCount, double[] values, int[] rowIndices, int[] columnPointers)
            : base(rowCount, columnCount, values, rowIndices, columnPointers)
        {
        }

        #region Public functions

        /// <inheritdoc />
        public override int DropZeros(double tolerance = 0.0)
        {
            Func<int, int, double, bool> func;

            if (tolerance <= 0.0)
            {
                func = (i, j, aij) =>
                {
                    return (aij != 0.0);
                };
            }
            else
            {
                func = (i, j, aij) =>
                {
                    return (Math.Abs(aij) > tolerance);
                };
            }

            return Keep(func);
        }

        /// <inheritdoc />
        public override int Keep(Func<int, int, double, bool> func)
        {
            int i, j, nz = 0;

            for (j = 0; j < columns; j++)
            {
                i = ColumnPointers[j];

                // Record new location of col j.
                ColumnPointers[j] = nz;

                for (; i < ColumnPointers[j + 1]; i++)
                {
                    if (func(RowIndices[i], j, Values[i]))
                    {
                        // Keep A(i,j).
                        Values[nz] = Values[i];
                        RowIndices[nz] = RowIndices[i];
                        nz++;
                    }
                }
            }

            // Record new nonzero count.
            ColumnPointers[columns] = nz;

            if (AutoTrimStorage)
            {
                // Remove extra space.
                this.Resize(0);
            }

            return nz;
        }

        /// <inheritdoc />
        public override double L1Norm()
        {
            int nz = this.NonZerosCount;

            double sum, norm = 0.0;

            for (int j = 0; j < columns; j++)
            {
                sum = 0.0;
                for (int i = ColumnPointers[j]; i < ColumnPointers[j + 1]; i++)
                {
                    sum += Math.Abs(Values[i]);
                }
                norm = Math.Max(norm, sum);
            }

            return norm;
        }
        
        /// <inheritdoc />
        public override double InfinityNorm()
        {
            int nz = this.NonZerosCount;

            double norm = 0.0;

            var work = new double[rows];

            for (int j = 0; j < columns; j++)
            {
                for (int i = ColumnPointers[j]; i < ColumnPointers[j + 1]; i++)
                {
                    work[RowIndices[i]] += Math.Abs(Values[i]);
                }
            }

            for (int j = 0; j < rows; j++)
            {
                norm = Math.Max(norm, work[j]);
            }

            return norm;
        }

        /// <inheritdoc />
        public override double FrobeniusNorm()
        {
            int nz = this.NonZerosCount;

            double sum = 0.0, norm = 0.0;
            
            for (int i = 0; i < nz; i++)
            {
                sum = Math.Abs(Values[i]);
                norm += sum * sum;
            }

            return Math.Sqrt(norm);
        }

        #endregion

        #region Linear Algebra (Vector)

        /// <inheritdoc />
        public override void Multiply(double[] x, double[] y)
        {
            var ax = this.Values;
            var ap = this.ColumnPointers;
            var ai = this.RowIndices;

            // Clear y.
            for (int i = 0; i < rows; i++)
            {
                y[i] = 0.0;
            }

            int end;

            for (int j = 0; j < columns; j++)
            {
                end = ap[j + 1];

                // Loop over the rows.
                for (int k = ap[j]; k < end; k++)
                {
                    y[ai[k]] += x[j] * ax[k];
                }
            }
        }

        /// <inheritdoc />
        public override void Multiply(double alpha, double[] x, double beta, double[] y)
        {
            var ax = this.Values;
            var ap = this.ColumnPointers;
            var ai = this.RowIndices;

            // Scale y by beta
            for (int j = 0; j < rows; j++)
            {
                y[j] = beta * y[j];
            }

            int end;
            double xi;

            for (int i = 0; i < columns; i++)
            {
                xi = alpha * x[i];

                end = ap[i + 1];

                for (int k = ap[i]; k < end; k++)
                {
                    y[ai[k]] += ax[k] * xi;
                }
            }
        }

        /// <inheritdoc />
        public override void TransposeMultiply(double[] x, double[] y)
        {
            var ax = this.Values;
            var ap = this.ColumnPointers;
            var ai = this.RowIndices;

            double yi;

            for (int i = 0; i < columns; i++)
            {
                yi = 0.0;

                // Compute the inner product of row i with vector x
                for (int k = ap[i]; k < ap[i + 1]; k++)
                {
                    yi += ax[k] * x[ai[k]];
                }

                // Store result in y(i) 
                y[i] = yi;
            }
        }

        /// <inheritdoc />
        public override void TransposeMultiply(double alpha, double[] x, double beta, double[] y)
        {
            var ax = this.Values;
            var ap = this.ColumnPointers;
            var ai = this.RowIndices;

            double yi;

            int end, start = ap[0];

            for (int i = 0; i < columns; i++)
            {
                end = ap[i + 1];

                yi = beta * y[i];
                for (int k = start; k < end; k++)
                {
                    yi += alpha * ax[k] * x[ai[k]];
                }
                y[i] = yi;

                start = end;
            }
        }

        #endregion

        #region Linear Algebra (Matrix)

        /// <inheritdoc />
        public override void Add(double alpha, double beta, CompressedColumnStorage<double> other,
            CompressedColumnStorage<double> result)
        {
            if (other == null)
            {
                throw new ArgumentNullException(nameof(other));
            }

            if (result == null)
            {
                throw new ArgumentNullException(nameof(result));
            }

            int p, j, nz = 0;

            int m = this.rows;
            int n = this.columns;

            // check inputs
            if (m != other.RowCount || n != other.ColumnCount)
            {
                throw new ArgumentException(Resources.MatrixDimensions, nameof(other));
            }

            // Workspace
            var w = new int[m];
            var x = new double[m];

            // Allocate result: (anz + bnz) is an upper bound

            var ci = result.ColumnPointers;
            var cj = result.RowIndices;
            var cx = result.Values;

            for (j = 0; j < n; j++)
            {
                ci[j] = nz; // column j of C starts here
                nz = this.Scatter(j, alpha, w, x, j + 1, result, nz); // alpha*A(:,j)
                nz = other.Scatter(j, beta, w, x, j + 1, result, nz); // beta*B(:,j)

                for (p = ci[j]; p < nz; p++)
                {
                    cx[p] = x[cj[p]];
                }
            }

            // Finalize the last column
            ci[n] = nz;

            if (AutoTrimStorage)
            {
                // Remove extra space.
                result.Resize(0);
            }

            Helper.SortIndices(result);
        }

        /// <inheritdoc />
        public override void Multiply(CompressedColumnStorage<double> other, CompressedColumnStorage<double> result)
        {
            if (other == null)
            {
                throw new ArgumentNullException(nameof(other));
            }

            if (result == null)
            {
                throw new ArgumentNullException(nameof(result));
            }

            int p, j, nz = 0;
            int[] cp, ci;
            double[] cx;

            int m = this.rows;
            int n = other.ColumnCount;

            int anz = this.NonZerosCount;
            int bnz = other.NonZerosCount;

            if (this.ColumnCount != other.RowCount)
            {
                throw new ArgumentException(Resources.MatrixDimensions, nameof(other));
            }

            if ((m > 0 && this.ColumnCount == 0) || (other.RowCount == 0 && n > 0))
            {
                throw new Exception(Resources.InvalidDimensions);
            }

            if (result.RowCount != m || result.ColumnCount != n)
            {
                throw new ArgumentException(Resources.InvalidDimensions, nameof(result));
            }

            var bp = other.ColumnPointers;
            var bi = other.RowIndices;
            var bx = other.Values;

            // Workspace
            var w = new int[m];
            var x = new double[m];

            cp = result.ColumnPointers;
            for (j = 0; j < n; j++)
            {
                if (nz + m > result.Values.Length)
                {
                    // Might throw out of memory exception.
                    result.Resize(2 * (result.Values.Length) + m);
                }
                ci = result.RowIndices;
                cx = result.Values; // C.i and C.x may be reallocated
                cp[j] = nz; // column j of C starts here
                for (p = bp[j]; p < bp[j + 1]; p++)
                {
                    nz = this.Scatter(bi[p], bx[p], w, x, j + 1, result, nz);
                }

                for (p = cp[j]; p < nz; p++)
                {
                    cx[p] = x[ci[p]];
                }
            }
            cp[n] = nz; // finalize the last column of C

            if (AutoTrimStorage)
            {
                // Remove extra space.
                result.Resize(0);
            }

            Helper.SortIndices(result);
        }

        /// <inheritdoc />
        public override CompressedColumnStorage<double> ParallelMultiply(CompressedColumnStorage<double> other, ParallelOptions options = null)
        {
            // Check inputs
            if (other == null)
            {
                throw new ArgumentNullException(nameof(other));
            }

            if (this.ColumnCount != other.RowCount)
            {
                throw new ArgumentException(Resources.MatrixDimensions);
            }

            int m = this.rows;
            int n = other.ColumnCount;

            if ((m > 0 && this.ColumnCount == 0) || (other.RowCount == 0 && n > 0))
            {
                throw new Exception(Resources.InvalidDimensions);
            }

            int processorCount = Environment.ProcessorCount;

            // Allow for at least 2 threads with 4 columns each
            if (m <= 0 || n < 2 * 4 || processorCount < 2)
            {
                return Multiply(other);
            }

            int anz = this.NonZerosCount;
            int bnz = other.NonZerosCount;

            // Heuristics to determine whether parallel multiplication is faster
            // Number of ops to exceed parallel overhead
            const int min_total_ops = 150000;

            // Total number of "x[i] += beta * Values[p]"
            long total_ops = (long)anz * bnz / this.ColumnCount;
            if (total_ops < min_total_ops)
            {
                return Multiply(other);
            }

            if (options == null)
            {
                options = new ParallelOptions() { MaxDegreeOfParallelism = processorCount };
            }
            else if (options.MaxDegreeOfParallelism < 0 || options.MaxDegreeOfParallelism > processorCount)
            {
                options.MaxDegreeOfParallelism = processorCount;
            }

            // With such a large overall threshold, there is no need of a per-thread threshold around 3000 ops.
            var nblocks = Math.Min(options.MaxDegreeOfParallelism, n);

            var bp = other.ColumnPointers;
            var bi = other.RowIndices;
            var bx = other.Values;

            var results = new SparseMatrix[nblocks];
            var indices = new int[nblocks];
            var nresults = 0;
            for (var j = 0; j < nblocks; j++)
            {
                var start = j * n / nblocks;
                var end = (j + 1) * n / nblocks;
                var bnz2 = bp[end] - bp[start];
                if (bnz2 != 0)
                {
                    indices[nresults] = start;
                    results[nresults++] = new SparseMatrix(m, end - start, anz + bnz2);
                }
            }
            Parallel.For(0, nresults, options,
                index =>
                {
                    var result = results[index];
                    var rnz = 0;

                    // Workspace
                    var w = new int[m];
                    var x = new double[m];

                    var rcp = result.ColumnPointers;
                    var nc = result.ColumnCount;
                    for (var j = 0; j < nc; j++)
                    {
                        if (rnz + m > result.Values.Length)
                        {
                            // Might throw out of memory exception.
                            result.Resize(2 * (result.Values.Length) + m);
                        }
                        var ci = result.RowIndices;
                        var cx = result.Values; // C.i and C.x may be reallocated
                        rcp[j] = rnz; // column j of C starts here
                        var j2 = j + indices[index];
                        for (var p = bp[j2]; p < bp[j2 + 1]; p++)
                        {
                            rnz = this.Scatter(bi[p], bx[p], w, x, j + 1, result, rnz);
                        }

                        for (var p = rcp[j]; p < rnz; p++)
                        {
                            cx[p] = x[ci[p]];
                        }
                    }

                    rcp[nc] = rnz; // finalize the last column of C
                    Helper.SortIndices(result);
                });

            int nz = 0;
            for (var j = 0; j < nresults; j++)
            {
                nz += results[j].NonZerosCount;
            }
            var values = new double[nz];
            var ri = new int[nz];
            var cp = new int[n + 1];
            nz = 0;
            var prev = 0;
            for (var j = 0; j < nresults; j++)
            {
                var start = indices[j];
                for (var k = prev; k < start; k++)
                {
                    cp[k] = nz;
                }
                var result = results[j];
                var rcp = result.ColumnPointers;
                var nc = result.ColumnCount;
                prev = start + nc;
                for (var k = 0; k < nc; k++)
                {
                    cp[start + k] = nz + rcp[k];
                }
                var rnz = result.NonZerosCount;
                Array.Copy(result.Values, 0, values, nz, rnz);
                Array.Copy(result.RowIndices, 0, ri, nz, rnz);
                nz += rnz;
            }
            for (var k = prev; k <= n; k++)
            {
                cp[k] = nz;
            }
            return new SparseMatrix(m, n, values, ri, cp);
        }

        #endregion

        /// <inheritdoc />
        public override bool Equals(Matrix<double> other, double tolerance)
        {
            var o = other as SparseMatrix;

            if (o == null)
            {
                return false;
            }

            int nz = this.NonZerosCount;

            if (this.columns != o.ColumnCount || this.rows != o.RowCount || nz != o.NonZerosCount)
            {
                return false;
            }

            for (int i = 0; i < this.columns; i++)
            {
                if (this.ColumnPointers[i] != o.ColumnPointers[i])
                {
                    return false;
                }
            }

            for (int i = 0; i < nz; i++)
            {
                if (this.RowIndices[i] != o.RowIndices[i])
                {
                    return false;
                }

                // TODO: should compare relative values!
                if (Math.Abs(this.Values[i] - o.Values[i]) > tolerance)
                {
                    return false;
                }
            }

            return true;
        }

        #region Internal methods
        
        internal override void Cleanup()
        {
            int i, j, p, q, nnz = 0;
            int[] marker = new int[rows];

            for (j = 0; j < rows; j++)
            {
                marker[j] = -1; // Row j not yet seen.
            }

            for (i = 0; i < columns; i++)
            {
                q = nnz; // Column i will start at q
                for (p = ColumnPointers[i]; p < ColumnPointers[i + 1]; p++)
                {
                    j = RowIndices[p]; // A(i,j) is nonzero
                    if (marker[j] >= q)
                    {
                        Values[marker[j]] += Values[p]; // A(i,j) is a duplicate
                    }
                    else
                    {
                        marker[j] = nnz; // Record where column j occurs
                        RowIndices[nnz] = j; // Keep A(i,j)
                        Values[nnz] = Values[p];

                        nnz += 1;
                    }
                }
                ColumnPointers[i] = q; // Record start of row i
            }

            this.ColumnPointers[columns] = nnz;

            // Remove extra space from arrays
            this.Resize(0);
        }

        /// <summary>
        /// Scatters and sums a sparse vector A(:,j) into a dense vector, x = x + beta * A(:,j).
        /// </summary>
        /// <param name="j">the column of A to use</param>
        /// <param name="beta">scalar multiplied by A(:,j)</param>
        /// <param name="w">size m, node i is marked if w[i] = mark</param>
        /// <param name="x">size m, not null</param>
        /// <param name="mark">mark value of w</param>
        /// <param name="mat">pattern of x accumulated in C.i</param>
        /// <param name="nz">pattern of x placed in C starting at C.i[nz]</param>
        /// <returns>new value of nz, -1 on error</returns>
        internal override int Scatter(int j, double beta, int[] w, double[] x, int mark,
            CompressedColumnStorage<double> mat, int nz)
        {
            int i, p;

            if (w == null || mat == null) return -1; // check inputs

            if (x == null)
            {
                throw new ArgumentNullException(nameof(x));
            }

            var cj = mat.RowIndices;

            for (p = ColumnPointers[j]; p < ColumnPointers[j + 1]; p++)
            {
                i = RowIndices[p]; // A(i,j) is nonzero

                if (w[i] < mark)
                {
                    w[i] = mark; // i is new entry in column j
                    x[i] = beta * Values[p]; // x(i) = beta*A(i,j)
                    cj[nz++] = i; // add i to pattern of C(:,j)
                }
                else
                {
                    x[i] += beta * Values[p]; // i exists in C(:,j) already
                }
            }

            return nz;
        }

        #endregion
    }
}
