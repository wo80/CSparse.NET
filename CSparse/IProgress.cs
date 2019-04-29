
namespace CSparse
{
    // If the project gets updated to .NET 4.5, this interface can be replaced
    // with the official System.IProgress<T> interface.

    /// <summary>
    /// Used to report progress of a factorization.
    /// </summary>
    public interface IProgress<T>
    {
        /// <summary>
        /// Reports a progress update.
        /// </summary>
        /// <param name="value">The value of the updated progress.</param>
        /// <remarks>
        /// For typeof(T) == double (as used in the factorizations), the value will range from 0.0 to 1.0.
        /// </remarks>
        void Report(T value);
    }
}
