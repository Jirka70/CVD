
using CVD.shape;

namespace CVD.util
{
    internal class Point3DGenerator
    {
        private static readonly Random RANDOM = new Random();
        private static readonly int Z_MAX_BOUND = 1000;

        private static Point3D GenerateRandomPoint(int windowWidth, int windowHeight)
        {
            int x = RANDOM.Next(windowWidth);
            int y = RANDOM.Next(windowHeight);
            int z = RANDOM.Next(Z_MAX_BOUND);

            return new Point3D(x, y, z);
        }

        public static List<Point3D> GenerateRandomPoints(int numberOfPoints, int windowWidth,
            int windowHeight)
        {
            List<Point3D> points = new();
            for (int i = 0; i < numberOfPoints; i++)
            {
                points.Add(GenerateRandomPoint(windowWidth, windowHeight));
            }
            return points;
        }
    }
}
