using VoronoiDiagrams;

namespace CVD
{
    internal class DelaunayTriangle
    {
        public readonly Point3D point1;
        public readonly Point3D point2;
        public readonly Point3D point3;

        public readonly Point3D circumCircleCenter;
        public readonly double circumCircleRadius;

        public DelaunayTriangle(Point3D point1, Point3D point2, Point3D point3)
        {
            Point3D[] pointArray = { point1, point2, point3 };
            Array.Sort(pointArray, new PointComparator());
            this.point1 = pointArray[0];
            this.point2 = pointArray[1];
            this.point3 = pointArray[2];
            this.circumCircleCenter = CalculateCircumCircleCenter();
            this.circumCircleRadius = CalculateEuclideanDistance(circumCircleCenter, point1);
        } 
        
        private Point3D CalculateCircumCircleCenter()
        {
            /*double ax = point1.X;
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
                    + (Math.Pow(cx, 2) + Math.Pow(cy, 2)) * (bx - ax));*/

            double acLength = point3.CalculateEuclideanDistance(point1);
            double abLength = point2.CalculateEuclideanDistance(point1);
            Point3D ab = point2.Subtract(point1);
            Point3D ac = point3.Subtract(point1);
            Point3D abXac = ab.CalculateCrossProduct(ac);


            Point3D point = abXac.Multiply(acLength*acLength).CalculateCrossProduct(ab).Add(ac.CalculateCrossProduct(abXac).Multiply(abLength*abLength));
            double denominator = 2f * Math.Pow(abXac.CalculateDistanceFromOriginOfCoordinateSystem(), 2);

   
            return point1.Add(new(point.X / denominator, point.Y / denominator, point.Z / denominator));
        }

        private double CalculateEuclideanDistance(Point3D a, Point3D b)
        {
            return Math.Sqrt(Math.Pow(a.X - b.X, 2) + Math.Pow(a.Y - b.Y, 2) + Math.Pow(a.Z - b.Z, 2));
        }

        public bool HasVertex(Point3D point)
        {
            return point1.Equals(point) || point2.Equals(point) || point3.Equals(point);
        }

        public bool HasEdge(Edge edge)
        {
            return HasVertex(edge.startingPoint) && HasVertex(edge.endingPoint);
        }

        public bool IsInCircumcircle(double x, double y, double z)
        {
            Point3D point = new(x, y, z);
            double distanceFromCircumcirleCenterToPoint = CalculateEuclideanDistance(circumCircleCenter, point);
            return distanceFromCircumcirleCenterToPoint < circumCircleRadius;
        }

        public override bool Equals(object? obj)
        {
            if (obj == null || GetType() != obj.GetType())
            {
                return false;
            }

            DelaunayTriangle other = (DelaunayTriangle) obj;
            return point1.Equals(other.point1) && point2.Equals(other.point2) && point3.Equals(other.point3);
        }

        public override int GetHashCode()
        {
            return HashCode.Combine(point1.X, point1.Y, point2.X, point2.Y, point3.X, point3.Y);
        }
    }
}
