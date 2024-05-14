using CVD.voronoi;

namespace CVD.util
{
    internal class ClockwisePointComparator : IComparer<VoronoiPoint2D>
    {
        private readonly VoronoiPoint2D center;
        public ClockwisePointComparator(VoronoiPoint2D center)
        {
            this.center = center;
        }
        public int Compare(VoronoiPoint2D a, VoronoiPoint2D b)
        {
            double angleA = Math.Atan2(a.Y - center.Y, a.X - center.X);
            double angleB = Math.Atan2(b.Y - center.Y, b.X - center.X);
            return angleA.CompareTo(angleB);
        }
    }
}
