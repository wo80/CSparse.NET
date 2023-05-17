
namespace CSparse.Tests.Complex
{
    using CSparse.Complex;
    using CSparse.Storage;
    using NUnit.Framework;
    using NUnit.Framework.Internal;

    using Complex = System.Numerics.Complex;

    [DefaultFloatingPointTolerance(1e-12)]
    public class SolverHelperTest
    {
        private readonly CompressedColumnStorage<Complex> L;
        private readonly CompressedColumnStorage<Complex> U;

        public SolverHelperTest()
        {
            var A = ResourceLoader.Get<Complex>("hermitian-40-spd.mat");

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

            Helper.AssertEqual(L.ColumnCount, x, b);
        }

        [Test]
        public void TestSolveLowerTranspose()
        {
            var x = Vector.Create(L.ColumnCount, 1.0);
            var b = Helper.TransposeMultiply(L, x);

            SolverHelper.SolveLowerTranspose(L, b);

            Helper.AssertEqual(L.ColumnCount, x, b);
        }

        [Test]
        public void TestSolveUpper()
        {
            var x = Vector.Create(U.ColumnCount, 1.0);
            var b = Helper.Multiply(U, x);

            SolverHelper.SolveUpper(U, b);

            Helper.AssertEqual(L.ColumnCount, x, b);
        }

        [Test]
        public void TestSolveUpperTranspose()
        {
            var x = Vector.Create(U.ColumnCount, 1.0);
            var b = Helper.TransposeMultiply(U, x);

            SolverHelper.SolveUpperTranspose(U, b);

            Helper.AssertEqual(L.ColumnCount, x, b);
        }
    }
}
