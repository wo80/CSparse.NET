
namespace CSparse.Tests
{
    using NUnit.Framework;

    public class PermutationTest
    {
        [Test]
        public void TestCreate()
        {
            var expected = new[] { 0, 1, 2 };
            var actual = Permutation.Create(3, 0);

            CollectionAssert.AreEqual(expected, actual);

            expected = new[] { 2, 1, 0 };
            actual = Permutation.Create(3, -1);

            CollectionAssert.AreEqual(expected, actual);
        }

        [Test]
        public void TestApply()
        {
            var vector = new double[] { 1d, 2d, 3d };
            var p = new[] { 2, 0, 1 };

            var expected = new double[] { 3d, 1d, 2d };
            var actual = new double[3];
            
            Permutation.Apply(p, vector, actual, 3);

            CollectionAssert.AreEqual(expected, actual);
        }

        [Test]
        public void TestApplyInverse()
        {
            var vector = new double[] { 1d, 2d, 3d };
            var p = new[] { 2, 0, 1 };

            var expected = new double[] { 2d, 3d, 1d };
            var actual = new double[3];

            Permutation.ApplyInverse(p, vector, actual, 3);

            CollectionAssert.AreEqual(expected, actual);
        }

        [Test]
        public void TestInvert()
        {
            var p = new[] { 2, 0, 1 };

            var expected = new[] { 1, 2, 0 };
            var actual = Permutation.Invert(p);

            CollectionAssert.AreEqual(expected, actual);
        }

        [Test]
        public void TestIsValid()
        {
            var p_valid = new[] { 2, 0, 1 };

            Assert.IsTrue(Permutation.IsValid(p_valid));

            var p_invalid = new[] { 2, 1, 1 };

            Assert.IsFalse(Permutation.IsValid(p_invalid));
        }
    }
}