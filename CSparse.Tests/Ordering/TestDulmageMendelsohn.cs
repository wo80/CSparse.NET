
namespace CSparse.Tests.Ordering
{
    using CSparse.Ordering;
    using Microsoft.VisualStudio.TestTools.UnitTesting;
    using System;

    [TestClass]
    public class TestDulmageMendelsohn
    {
        [TestMethod]
        public void TestGenerate1()
        {
            // Load matrix from a file.
            var A = ResourceLoader.Get<double>("general-40x40.mat");

            int n = A.ColumnCount;

            // Randomized Dulmage-Mendelsohn analysis.
            var dm = DulmageMendelsohn.Generate(A, 1);

            Assert.IsTrue(dm.StructuralRank == n);
        }

        [TestMethod]
        public void TestGenerate2()
        {
            // Load matrix from a file.
            var A = ResourceLoader.Get<double>("general-40x20.mat");

            int n = Math.Min(A.RowCount, A.ColumnCount);

            // Randomized Dulmage-Mendelsohn analysis.
            var dm = DulmageMendelsohn.Generate(A, 1);

            Assert.IsTrue(dm.StructuralRank == n);
        }
    }
}
