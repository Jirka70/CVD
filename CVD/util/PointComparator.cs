using CVD.voronoi;

namespace CVD.util
{
    internal class PointComparator : IComparer<VoronoiPoint2D>
    {
        public int Compare(VoronoiPoint2D p1, VoronoiPoint2D p2)
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
