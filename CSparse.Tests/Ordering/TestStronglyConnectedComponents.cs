
namespace CSparse.Tests.Ordering
{
    using CSparse.Double;
    using CSparse.Ordering;
    using NUnit.Framework;
    using System;

    public class TestStronglyConnectedComponents
    {
        [Test]
        public void TestScc()
        {
            // https://en.wikipedia.org/wiki/Strongly_connected_component
            //
            // Adjacency list of directed graph with 8 nodes:
            //
            // [1] -> 2
            // [2] -> 3 -> 5 -> 6
            // [3] -> 4 -> 7
            // [4] -> 3 -> 8
            // [5] -> 1 -> 6
            // [6] -> 7
            // [7] -> 6
            // [8] -> 4 -> 7
            //
            // Adjacency matrix:
            
            var A = SparseMatrix.OfRowMajor(8, 8, new double[8 * 8]
            {
                1, 1, 0, 0, 0, 0, 0, 0,
                0, 1, 1, 0, 1, 1, 0, 0,
                0, 0, 1, 1, 0, 0, 1, 0,
                0, 0, 1, 1, 0, 0, 0, 1,
                1, 0, 0, 0, 1, 1, 0, 0,
                0, 0, 0, 0, 0, 1, 1, 0,
                0, 0, 0, 0, 0, 1, 1, 0,
                0, 0, 0, 1, 0, 0, 1, 1
            });

            int n = A.ColumnCount;
            
            var scc = StronglyConnectedComponents.Generate(A);

            var r = new int[] { 0, 3, 6, 8 };
            var p = new int[] { 0, 1, 4, 2, 3, 7, 5, 6 };

            Assert.AreEqual(scc.Blocks, 3);
            CollectionAssert.AreEqual(scc.BlockPointers, r);
            CollectionAssert.AreEqual(scc.Indices, p);
        }
    }
}
