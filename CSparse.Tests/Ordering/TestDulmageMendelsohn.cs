namespace CSparse.Tests.Ordering
{
    using CSparse.Double;
    using CSparse.Ordering;
    using CSparse.Storage;
    using NUnit.Framework;
    using System;

    public class TestDulmageMendelsohn
    {
        [Test]
        public void TestGenerate1()
        {
            // Load matrix from a file.
            var A = ResourceLoader.Get<double>("general-40x40.mat");

            int n = A.ColumnCount;

            // Randomized Dulmage-Mendelsohn analysis.
            var dm = DulmageMendelsohn.Generate(A, 1);

            Assert.That(dm.StructuralRank == n, Is.True);
        }

        [Test]
        public void TestGenerate2()
        {
            // Load matrix from a file.
            var A = ResourceLoader.Get<double>("general-40x20.mat");

            int n = Math.Min(A.RowCount, A.ColumnCount);

            // Randomized Dulmage-Mendelsohn analysis.
            var dm = DulmageMendelsohn.Generate(A, 1);

            Assert.That(dm.StructuralRank == n, Is.True);
        }

        [Test]
        public void TestGenerate3()
        {
            var A = SparseMatrix.OfRowMajor(8, 8,
            [
                11, 12,  0,  0,  0,  0,  0,  0,
                 0, 22, 23,  0, 25, 26,  0,  0,
                 0,  0, 33, 34,  0,  0, 37,  0,
                 0,  0, 43, 44,  0,  0,  0, 48,
                51,  0,  0,  0, 55, 56,  0,  0,
                 0,  0,  0,  0,  0, 66, 67,  0,
                 0,  0,  0,  0,  0, 76, 77,  0,
                 0,  0,  0, 84,  0,  0, 87, 88
            ]);

            var S = SymbolicColumnStorage.Create(A);

            var dm = DulmageMendelsohn.Generate(S, 1);

            Assert.That(dm.StructuralRank, Is.EqualTo(8));
            Assert.That(dm.Blocks, Is.EqualTo(3));

            int[] expected = [0, 1, 4, 2, 3, 7, 5, 6];

            Assert.That(dm.RowPermutation, Is.EqualTo(expected).AsCollection);
            Assert.That(dm.ColumnPermutation, Is.EqualTo(expected).AsCollection);

            expected = [0, 3, 6, 8];

            Assert.That(dm.BlockRowPointers, Is.EqualTo(expected).AsCollection);
            Assert.That(dm.BlockColumnPointers, Is.EqualTo(expected).AsCollection);

            expected = [0, 0, 8, 8, 8];

            Assert.That(dm.CoarseRowDecomposition, Is.EqualTo(expected).AsCollection);

            expected = [0, 0, 0, 8, 8];

            Assert.That(dm.CoarseColumnDecomposition, Is.EqualTo(expected).AsCollection);
        }
    }
}
