
namespace CSparse.Tests.Double
{
    using System;
    using System.Globalization;
    using System.IO;

    class DenseTestDataReader
    {
        public static DenseTestData<double> Read(Stream stream)
        {
            using (var reader = new StreamReader(stream))
            {
                string line, name, value;

                line = reader.ReadLine();

                GetItem(line, out name, out value);

                if (name != "size")
                {
                    throw new FormatException("Expected first line = size.");
                }

                var data = ReadSize(value);

                int m = data.RowCount;
                int n = data.ColumnCount;

                while ((line = reader.ReadLine()) != null)
                {
                    GetItem(line, out name, out value);

                    if (name == "A")
                    {
                        data.A = ReadMatrix(value, m, n);
                    }
                    else if (name == "B")
                    {
                        data.B = ReadMatrix(value, m, n);
                    }
                    else if (name == "x")
                    {
                        data.x = ReadVector(value, n);
                    }
                    else if (name == "y")
                    {
                        data.y = ReadVector(value, m);
                    }
                    else if (name == "A'")
                    {
                        data.AT = ReadMatrix(value, n, m);
                    }
                    else if (name == "B'")
                    {
                        data.BT = ReadMatrix(value, n, m);
                    }
                    else if (name == "A+B")
                    {
                        data.ApB = ReadMatrix(value, m, n);
                    }
                    else if (name == "A*B'")
                    {
                        data.AmBT = ReadMatrix(value, m, m);
                    }
                    else if (name == "A'*B")
                    {
                        data.ATmB = ReadMatrix(value, n, n);
                    }
                    else if (name == "A*x")
                    {
                        data.Ax = ReadVector(value, m);
                    }
                    else if (name == "A'*y")
                    {
                        data.ATy = ReadVector(value, n);
                    }
                    else if (name == "x'*B'")
                    {
                        data.xTBT = ReadVector(value, m);
                    }
                }

                return data;
            }
        }

        private static void GetItem(string line, out string name, out string value)
        {
            if (string.IsNullOrEmpty(line))
            {
                name = value = null;

                return;
            }

            int i = line.IndexOf(':');

            if (i < 0)
            {
                throw new FormatException("Expected name-value pair separated by colon.");
            }

            name = line.Substring(0, i).Trim();
            value = line.Substring(i + 1).Trim();
        }

        private static DenseTestData<double> ReadSize(string line)
        {
            var size = line.Split();

            return new DenseTestData<double>()
            {
                RowCount = int.Parse(size[0]),
                ColumnCount = int.Parse(size[1])
            };
        }

        private static double[,] ReadMatrix(string line, int m, int n)
        {
            var rows = line.Trim().Split(';');

            if (rows.Length != m)
            {
                throw new InvalidDataException();
            }

            var result = new double[m, n];

            for (int i = 0; i < m; i++)
            {
                var row = rows[i].Trim().Split();

                if (row.Length != n)
                {
                    throw new InvalidDataException();
                }

                for (int j = 0; j < n; j++)
                {
                    result[i, j] = double.Parse(row[j], CultureInfo.InvariantCulture);
                }
            }

            return result;
        }

        private static double[] ReadVector(string line, int size)
        {
            line = line.Replace(";", string.Empty).Trim();

            var values = line.Split(' ');

            if (values.Length != size)
            {
                throw new InvalidDataException();
            }

            var result = new double[size];

            for (int i = 0; i < size; i++)
            {
                result[i] = double.Parse(values[i], CultureInfo.InvariantCulture);
            }

            return result;
        }
    }
}
