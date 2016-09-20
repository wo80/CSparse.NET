
namespace CSparse
{
    // If the project gets updated to .NET 4.5, this interface can be replaced
    // with the official System.IProgress<T> interface.

    /// <summary>
    /// Used to report progress of a factorization.
    /// </summary>
    public interface IProgress
    {
        /// <summary>
        /// Reports a progress update.
        /// </summary>
        /// <param name="value">The value of the updated progress (from 0.0 to 1.0).</param>
        void Report(double value);
    }
}
