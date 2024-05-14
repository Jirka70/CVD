namespace CVD.voronoi
{
    internal class VoronoiPoint2D
    {

        public readonly double X;
        public readonly double Y;

        public VoronoiPoint2D(double x, double y)
        {
            X = x;
            Y = y;
        }

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
            return HashCode.Combine(X, Y);
        }
    }
}
