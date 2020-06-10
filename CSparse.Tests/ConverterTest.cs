using CSparse.Storage;
using NUnit.Framework;

namespace CSparse.Tests
{
    public class ConverterTest
    {
        double[][] Jarray;
        double[,] Marray;
        CompressedColumnStorage<double> Jcsc = null;
        CompressedColumnStorage<double> Mcsc = null;

        [SetUp]
        public void Initialize()
        {
            // jagged array
            Jarray = new double[3][];
            for (int i = 0; i < 3; i++)
                Jarray[i] = new double[3];

            Jarray[0][1] = 2.0; Jarray[1][0] = 2.0;   // | 0 2 4 |
            Jarray[0][2] = 4.0; Jarray[2][0] = 4.0;   // | 2 0 0 |
            Jarray[2][2] = 6.0;                       // | 4 0 6 |
            Jcsc = CompressedColumnStorage<double>.OfJaggedArray(Jarray);

            // multidimentional array
            Marray = new double[,] 
            {
                {0, 2, 4 },
                {2, 0, 0 },
                {4, 0, 6 }
            };
            Mcsc = CompressedColumnStorage<double>.OfArray(Marray);
        }

        [Test]
        public void TestConversion()
        {
            Assert.NotNull(Jcsc);
            Assert.NotNull(Mcsc);

            int[] colPointers = new[] { 0, 2, 3, 5 };
            int[] rowIndices = new[] { 1, 2, 0, 0, 2 };
            double[] vals = new[] { 2.0, 4.0, 2.0, 4.0, 6.0 };

            // jagged
            CollectionAssert.AreEqual(colPointers, Jcsc.ColumnPointers);
            CollectionAssert.AreEqual(rowIndices, Jcsc.RowIndices);
            CollectionAssert.AreEqual(vals, Jcsc.Values);

            // multi-
            CollectionAssert.AreEqual(colPointers, Mcsc.ColumnPointers);
            CollectionAssert.AreEqual(rowIndices, Mcsc.RowIndices);
            CollectionAssert.AreEqual(vals, Mcsc.Values);
        }

        [Test]
        public void TestEquality()
        {
            // both should be equal after conversion
            Assert.AreEqual(Jcsc, Mcsc);
        }
    }
}
