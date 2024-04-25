namespace VoronoiDiagrams
{
    internal class VoronoiPoint2D(double x, double y)
    {
        public readonly double X = x;
        public readonly double Y = y;

        public override bool Equals(object? obj)
        {
            if (obj == null || GetType() != obj.GetType())
            {
                return false;
            }

            VoronoiPoint2D other = (VoronoiPoint2D)obj;
            return X == other.X && Y == other.Y;
        }

        public override int GetHashCode()
        {
            return HashCode.Combine(x, y);
        }
    }
}
