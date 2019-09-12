
namespace CSparse.Tests.Complex
{
    using System;
    using System.Collections;
    using System.Collections.Generic;
    using System.Numerics;

    // Setting DefaultFloatingPointTolerance to 1e-8 seems to work only for double/float types.
    // For System.Numerics.Complex NUnit uses the default "Equals" implementation.

    public class ComplexNumberComparer : IComparer, IComparer<Complex>
    {
        public static ComplexNumberComparer Default = new ComplexNumberComparer();

        // Floating point tolerance
        const double TOL = 1e-8;

        public int Compare(object x, object y)
        {
            return Compare((Complex)x, (Complex)y);
        }

        public int Compare(Complex x, Complex y)
        {
            return (Math.Abs(x.Real - y.Real) < TOL && Math.Abs(x.Imaginary - y.Imaginary) < TOL) ? 0 : 1;
        }
    }
}
