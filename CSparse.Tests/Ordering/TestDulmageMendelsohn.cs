namespace CSparse.Tests.Ordering
{
    using CSparse.Ordering;
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
    }
}
