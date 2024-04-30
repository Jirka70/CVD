using VoronoiDiagrams;

namespace CVD
{
    internal class PointComparator : IComparer<Point3D>
    {
        public int Compare(Point3D p1, Point3D p2)
        {
            int xCordComparison = p1.X.CompareTo(p2.X);
            if (xCordComparison != 0)
            {
                return xCordComparison;
            }

            return p1.Y.CompareTo(p2.Y);
        }
    }
}
