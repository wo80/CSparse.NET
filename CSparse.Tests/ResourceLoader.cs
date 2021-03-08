using CSparse.IO;
using CSparse.Storage;
using System;
using System.Collections.Generic;
using System.IO;
using System.Reflection;

namespace CSparse.Tests
{
    static class ResourceLoader
    {
        private const string NS = "CSparse.Tests.{Type}.Data";

        private static Dictionary<string, object> cache = new Dictionary<string, object>();

        public static Stream GetStream(string resource, string type)
        {
            string path = NS.Replace("{Type}", type) + "." + resource;

            return Assembly.GetExecutingAssembly().GetManifestResourceStream(path);
        }

        public static CompressedColumnStorage<T> Get<T>(string resource)
            where T : struct, IEquatable<T>, IFormattable
        {
            try
            {
                string type = (typeof(T) == typeof(double)) ? "Double" : "Complex";

                string path = NS.Replace("{Type}", type) + "." + resource;

                object obj;

                if (cache.TryGetValue(path, out obj))
                {
                    return (CompressedColumnStorage<T>)obj;
                }

                var stream = Assembly.GetExecutingAssembly().GetManifestResourceStream(path);

                var matrix = MatrixMarketReader.ReadMatrix<T>(stream);

                cache.Add(path, matrix);

                return matrix;
            }
            catch
            {
                throw;
            }
        }
    }
}
