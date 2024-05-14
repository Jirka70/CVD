using CVD.shape;
using CVD.util;
using CVD.voronoi;

namespace CVD.delaunay
{
    internal class DelaunayTriangle
    {
        public readonly VoronoiPoint2D point1;
        public readonly VoronoiPoint2D point2;
        public readonly VoronoiPoint2D point3;

        public readonly VoronoiPoint2D circumCircleCenter;
        public readonly double circumCircleRadius;

        public DelaunayTriangle(VoronoiPoint2D point1, VoronoiPoint2D point2, VoronoiPoint2D point3)
        {
            VoronoiPoint2D[] pointArray = { point1, point2, point3 };
            Array.Sort(pointArray, new PointComparator());
            this.point1 = pointArray[0];
            this.point2 = pointArray[1];
            this.point3 = pointArray[2];
            circumCircleCenter = CalculateCircumCircleCenter();
            circumCircleRadius = CalculateEuclideanDistance(circumCircleCenter, point1);
        }

        private VoronoiPoint2D CalculateCircumCircleCenter()
        {
            double ax = point1.X;
            double bx = point2.X;
            double cx = point3.X;
            double ay = point1.Y;
            double by = point2.Y;
            double cy = point3.Y;
            double d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));

            double centerX = 1 / d * ((Math.Pow(ax, 2) + Math.Pow(ay, 2)) * (by - cy) + (Math.Pow(bx, 2)
                + Math.Pow(by, 2)) * (cy - ay)
                + (Math.Pow(cx, 2) + Math.Pow(cy, 2)) * (ay - by));

            double centerY = 1 / d * ((Math.Pow(ax, 2) + Math.Pow(ay, 2)) * (cx - bx)
                    + (Math.Pow(bx, 2) + Math.Pow(by, 2)) * (ax - cx)
                    + (Math.Pow(cx, 2) + Math.Pow(cy, 2)) * (bx - ax));

            return new(centerX, centerY);
        }

        private double CalculateEuclideanDistance(VoronoiPoint2D a, VoronoiPoint2D b)
        {
            return Math.Sqrt(Math.Pow(a.X - b.X, 2) + Math.Pow(a.Y - b.Y, 2));
        }

        public bool HasVertex(VoronoiPoint2D point)
        {
            return point1.Equals(point) || point2.Equals(point) || point3.Equals(point);
        }

        public bool HasEdge(Edge edge)
        {
            return HasVertex(edge.startingPoint) && HasVertex(edge.endingPoint);
        }

        public bool IsInCircumcircle(double x, double y)
        {
            VoronoiPoint2D point = new(x, y);
            double distanceFromCircumcirleCenterToPoint = CalculateEuclideanDistance(circumCircleCenter, point);
            return distanceFromCircumcirleCenterToPoint < circumCircleRadius;
        }

        public override bool Equals(object? obj)
        {
            if (obj == null || GetType() != obj.GetType())
            {
                return false;
            }

            DelaunayTriangle other = (DelaunayTriangle)obj;
            return point1.Equals(other.point1) && point2.Equals(other.point2) && point3.Equals(other.point3);
        }

        public override int GetHashCode()
        {
            return HashCode.Combine(point1.X, point1.Y, point2.X, point2.Y, point3.X, point3.Y);
        }
    }
}
