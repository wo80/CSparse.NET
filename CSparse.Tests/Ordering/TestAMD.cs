using CSparse.Ordering;
using CSparse.Storage;
using NUnit.Framework;

namespace CSparse.Tests.Ordering
{
    class TestAMD
    {
        [Test]
        public void TestAMD1()
        {
            int[] ap = [0, 9, 15, 21, 27, 33, 39, 48, 57, 61, 70, 76, 82, 88, 94, 100, 106, 110, 119, 128, 137, 143, 152, 156, 160];
            int[] ai = [
            /* col  0 */  0, 5, 6, 12, 13, 17, 18, 19, 21,
            /* col  1 */  1, 8, 9, 13, 14, 17,
            /* col  2 */  2, 6, 11, 20, 21, 22,
            /* col  3 */  3, 7, 10, 15, 18, 19,
            /* col  4 */  4, 7, 9, 14, 15, 16,
            /* col  5 */  0, 5, 6, 12, 13, 17,
            /* col  6 */  0, 2, 5, 6, 11, 12, 19, 21, 23,
            /* col  7 */  3, 4, 7, 9, 14, 15, 16, 17, 18,
            /* col  8 */  1, 8, 9, 14,
            /* col  9 */  1, 4, 7, 8, 9, 13, 14, 17, 18,
            /* col 10 */  3, 10, 18, 19, 20, 21,
            /* col 11 */  2, 6, 11, 12, 21, 23,
            /* col 12 */  0, 5, 6, 11, 12, 23,
            /* col 13 */  0, 1, 5, 9, 13, 17,
            /* col 14 */  1, 4, 7, 8, 9, 14,
            /* col 15 */  3, 4, 7, 15, 16, 18,
            /* col 16 */  4, 7, 15, 16,
            /* col 17 */  0, 1, 5, 7, 9, 13, 17, 18, 19,
            /* col 18 */  0, 3, 7, 9, 10, 15, 17, 18, 19,
            /* col 19 */  0, 3, 6, 10, 17, 18, 19, 20, 21,
            /* col 20 */  2, 10, 19, 20, 21, 22,
            /* col 21 */  0, 2, 6, 10, 11, 19, 20, 21, 22,
            /* col 22 */  2, 20, 21, 22,
            /* col 23 */  6, 11, 12, 23 ];

            var S = new SymbolicColumnStorage(24, 24, ap, ai, false);

            var p = AMD.Generate(S, ColumnOrdering.MinimumDegreeAtA);

            int[] expected = [8, 16, 4, 14, 15, 1, 7, 9, 22, 23, 2, 11, 3, 12, 13, 17, 18, 0, 10, 19, 20, 6, 21, 5, 24];

            Assert.That(Permutation.IsValid(p), Is.True);
            Assert.That(p, Is.EqualTo(expected).AsCollection);
        }

        [Test]
        public void TestAMD2()
        {
            int[] ap = [0, 9, 14, 20, 28, 33, 37, 44, 53, 58, 63, 63, 66, 69, 72, 75, 78, 82, 86, 91, 97, 101, 112, 112, 116];
            int[] ai = [
            /* col  0 */  0, 17, 18, 21, 5, 12, 5, 0, 13,
            /* col  1 */  14, 1, 8, 13, 17,
            /* col  2 */  2, 20, 11, 6, 11, 22,
            /* col  3 */  3, 3, 10, 7, 18, 18, 15, 19,
            /* col  4 */  7, 9, 15, 14, 16,
            /* col  5 */  5, 13, 6, 17,
            /* col  6 */  5, 0, 11, 6, 12, 6, 23,
            /* col  7 */  3, 4, 9, 7, 14, 16, 15, 17, 18,
            /* col  8 */  1, 9, 14, 14, 14,
            /* col  9 */  7, 13, 8, 1, 17,
            /* col 10 */
            /* col 11 */  2, 12, 23,
            /* col 12 */  5, 11, 12,
            /* col 13 */  0, 13, 17,
            /* col 14 */  1, 9, 14,
            /* col 15 */  3, 15, 16,
            /* col 16 */  16, 4, 4, 15,
            /* col 17 */  13, 17, 19, 17,
            /* col 18 */  15, 17, 19, 9, 10,
            /* col 19 */  17, 19, 20, 0, 6, 10,
            /* col 20 */  22, 10, 20, 21,
            /* col 21 */  6, 2, 10, 19, 20, 11, 21, 22, 22, 22, 22,
            /* col 22 */
            /* col 23 */  12, 11, 12, 23 ];

            var S = new SymbolicColumnStorage(24, 24, ap, ai, false);

            var p = AMD.Generate(S, ColumnOrdering.MinimumDegreeAtA);

            int[] expected = [10, 11, 23, 12, 2, 6, 8, 14, 15, 16, 4, 1, 9, 7, 18, 3, 5, 17, 0, 19, 20, 21, 13, 22, 24];

            Assert.That(Permutation.IsValid(p), Is.True);
            Assert.That(p, Is.EqualTo(expected).AsCollection);
        }
    }
}
