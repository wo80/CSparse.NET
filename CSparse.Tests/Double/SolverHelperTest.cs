
namespace CSparse.Tests.Double
{
    using CSparse.Double;
    using CSparse.Storage;
    using NUnit.Framework;
    using NUnit.Framework.Internal;

    [DefaultFloatingPointTolerance(1e-12)]
    public class SolverHelperTest
    {
        private readonly CompressedColumnStorage<double> L;
        private readonly CompressedColumnStorage<double> U;

        public SolverHelperTest()
        {
            var A = ResourceLoader.Get<double>("symmetric-40-spd.mat");

            L = A.Clone();
            L.Keep((i, j, a) => i >= j);

            U = A.Clone();
            U.Keep((i, j, a) => i <= j);
        }

        [Test]
        public void TestSolveLower()
        {
            var x = Vector.Create(L.ColumnCount, 1.0);
            var b = Helper.Multiply(L, x);

            SolverHelper.SolveLower(L, b);

            CollectionAssert.AreEqual(x, b);
        }

        [Test]
        public void TestSolveLowerTranspose()
        {
            var x = Vector.Create(L.ColumnCount, 1.0);
            var b = Helper.TransposeMultiply(L, x);

            SolverHelper.SolveLowerTranspose(L, b);

            CollectionAssert.AreEqual(x, b);
        }

        [Test]
        public void TestSolveUpper()
        {
            var x = Vector.Create(U.ColumnCount, 1.0);
            var b = Helper.Multiply(U, x);

            SolverHelper.SolveUpper(U, b);

            CollectionAssert.AreEqual(x, b);
        }

        [Test]
        public void TestSolveUpperTranspose()
        {
            var x = Vector.Create(U.ColumnCount, 1.0);
            var b = Helper.TransposeMultiply(U, x);

            SolverHelper.SolveUpperTranspose(U, b);

            CollectionAssert.AreEqual(x, b);
        }
    }
}
