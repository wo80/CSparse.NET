namespace CSparse.IO
{
    using CSparse.Storage;
    using System;
    using System.Globalization;
    using System.IO;
    using System.Numerics;

    internal enum MatrixMarketSymmetry
    {
        General,
        Symmetric,
        SkewSymmetric,
        Hermitian
    }

    /// <summary>
    /// Read files in Matrix Market format (only supports coordinate storage).
    /// </summary>
    public class MatrixMarketReader
    {
        private static readonly NumberFormatInfo nfi = CultureInfo.InvariantCulture.NumberFormat;
        private static readonly char[] separator = new char[] { ' ', '\t' };

        /// <summary>
        /// Read a matrix from file.
        /// </summary>
        public static CompressedColumnStorage<T> ReadMatrix<T>(string filePath)
            where T : struct, IEquatable<T>, IFormattable
        {
            using (var stream = File.OpenRead(filePath))
            {
                return ReadMatrix<T>(stream);
            }
        }

        /// <summary>
        /// Read a matrix from stream.
        /// </summary>
        public static CompressedColumnStorage<T> ReadMatrix<T>(Stream stream)
            where T : struct, IEquatable<T>, IFormattable
        {
            using (var reader = new StreamReader(stream))
            {
                return Converter.ToCompressedColumnStorage(ReadStorage<T>(reader));
            }
        }

        /// <summary>
        /// Read coordinate storage from a text reader.
        /// </summary>
        public static CoordinateStorage<T> ReadStorage<T>(TextReader reader, bool autoExpand = true)
            where T : struct, IEquatable<T>, IFormattable
        {
            bool complex, sparse, pattern;
            MatrixMarketSymmetry symmetry;
            ExpectHeader(reader, true, out complex, out sparse, out pattern, out symmetry);

            var parse = CreateValueParser<T>(complex, pattern);

            var sizes = ExpectLine(reader).Split(separator, StringSplitOptions.RemoveEmptyEntries);

            int rows = int.Parse(sizes[0]);
            int cols = int.Parse(sizes[1]);

            if (sparse)
            {
                int nnz = int.Parse(sizes[2]);

                var storage = new CoordinateStorage<T>(rows, cols, nnz);

                var s = BuildStorage(reader, storage, parse);

                bool lower = s.Item1;
                bool upper = s.Item2;

                // The Matrix Market specification states that symmetric/Hermitian matrices
                // should only store the lower triangular part, so we have to expand the
                // storage manually:
                if (autoExpand && symmetry != MatrixMarketSymmetry.General && !upper)
                {
                    ExpandStorage(symmetry, storage);
                }

                return storage;
            }

            throw new NotSupportedException();
        }

        private static Tuple<bool, bool> BuildStorage<T>(TextReader reader, CoordinateStorage<T> storage,
            Func<int, string[], T> parse)
            where T : struct, IEquatable<T>, IFormattable
        {
            int i, j;

            bool lower = false;
            bool upper = false;

            string line;
            while ((line = reader.ReadLine()) != null)
            {
                var trim = line.Trim();
                if (trim.Length > 0 && !trim.StartsWith("%"))
                {
                    var tokens = trim.Split(separator, StringSplitOptions.RemoveEmptyEntries);

                    i = int.Parse(tokens[0]) - 1;
                    j = int.Parse(tokens[1]) - 1;

                    storage.At(i, j, parse(2, tokens));

                    if (i > j)
                    {
                        lower = true;
                    }

                    if (i < j)
                    {
                        upper = true;
                    }
                }
            }

            return new Tuple<bool, bool>(lower, upper);
        }

        static void ExpandStorage<T>(MatrixMarketSymmetry symmetry, CoordinateStorage<T> storage)
            where T : struct, IEquatable<T>, IFormattable
        {
            var map = CreateSymmetryMap<T>(symmetry);

            var ai = storage.RowIndices;
            var aj = storage.ColumnIndices;
            var ax = storage.Values;

            int i, j, n = storage.NonZerosCount;

            for (int k = 0; k < n; k++)
            {
                i = ai[k];
                j = aj[k];

                if (i > j)
                {
                    storage.At(j, i, map(ax[k]));
                }
            }
        }

        #region Helper

        static void ExpectHeader(TextReader reader, bool matrix, out bool complex, out bool sparse, out bool pattern, out MatrixMarketSymmetry symmetry)
        {
            complex = pattern = false;
            string line;
            while ((line = reader.ReadLine()) != null)
            {
                line = line.Trim();
                if (line.StartsWith("%%MatrixMarket"))
                {
                    var tokens = line.ToLowerInvariant().Substring(15).Split(separator, StringSplitOptions.RemoveEmptyEntries);
                    if (tokens.Length < 2)
                    {
                        throw new FormatException(@"Expected MatrixMarket Header with 2-4 attributes: object format [field] [symmetry]; see http://math.nist.gov/MatrixMarket/ for details.");
                    }

                    if (tokens[0] != (matrix ? "matrix" : "vector"))
                    {
                        throw new FormatException("Expected matrix content.");
                    }

                    switch (tokens[1])
                    {
                        case "array":
                            sparse = false;
                            break;
                        case "coordinate":
                            sparse = true;
                            break;
                        default:
                            throw new NotSupportedException("Format type not supported.");
                    }

                    if (tokens.Length < 3)
                    {
                        complex = false;
                    }
                    else
                    {
                        switch (tokens[2])
                        {
                            case "real":
                            case "double":
                            case "integer":
                                complex = false;
                                break;
                            case "complex":
                                complex = true;
                                break;
                            case "pattern":
                                pattern = true;
                                break;
                            default:
                                throw new NotSupportedException("Field type not supported.");
                        }
                    }

                    if (tokens.Length < 4)
                    {
                        symmetry = MatrixMarketSymmetry.General;
                    }
                    else
                    {
                        switch (tokens[3])
                        {
                            case "general":
                                symmetry = MatrixMarketSymmetry.General;
                                break;
                            case "symmetric":
                                symmetry = MatrixMarketSymmetry.Symmetric;
                                break;
                            case "skew-symmetric":
                                symmetry = MatrixMarketSymmetry.SkewSymmetric;
                                break;
                            case "hermitian":
                                symmetry = MatrixMarketSymmetry.Hermitian;
                                break;
                            default:
                                throw new NotSupportedException("Symmetry type not supported");
                        }
                    }

                    return;
                }
            }

            throw new FormatException(@"Expected MatrixMarket Header, see http://math.nist.gov/MatrixMarket/ for details.");
        }

        static string ExpectLine(TextReader reader)
        {
            string line;
            while ((line = reader.ReadLine()) != null)
            {
                var trim = line.Trim();
                if (trim.Length > 0 && !trim.StartsWith("%"))
                {
                    return trim;
                }
            }

            throw new FormatException(@"End of file reached unexpectedly.");
        }

        static Func<int, string[], T> CreateValueParser<T>(bool sourceIsComplex, bool pattern)
        {
            if (pattern)
            {
                return (offset, tokens) => (T)(object)1.0;
            }

            if (typeof(T) == typeof(double))
            {
                // ignore imaginary part if source is complex
                return (offset, tokens) => (T)(object)double.Parse(tokens[offset], NumberStyles.Any, nfi);
            }

            if (typeof(T) == typeof(Complex))
            {
                return sourceIsComplex
                    ? ((offset, tokens) => (T)(object)new Complex(double.Parse(tokens[offset], NumberStyles.Any, nfi), double.Parse(tokens[offset + 1], NumberStyles.Any, nfi)))
                    : (Func<int, string[], T>)((offset, tokens) => (T)(object)new Complex(double.Parse(tokens[offset], NumberStyles.Any, nfi), 0d));
            }

            throw new NotSupportedException();
        }

        static Func<T, T> CreateSymmetryMap<T>(MatrixMarketSymmetry symmetry)
        {
            if (symmetry != MatrixMarketSymmetry.Hermitian)
            {
                return x => x;
            }

            if (typeof(T) == typeof(double))
            {
                return x => x;
            }

            if (typeof(T) == typeof(Complex))
            {
                return x => (T)(object)Complex.Conjugate((Complex)(object)x);
            }

            throw new NotSupportedException();
        }

        #endregion
    }
}
