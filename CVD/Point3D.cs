
using CVD;

namespace VoronoiDiagrams
{
   
    internal class Point3D {
        public readonly double X;
        public readonly double Y;
        public readonly double Z;

        public Point3D(double x, double y, double z)
        {
            X = x;
            Y = y;
            Z = z;
        }

        public double CalculateEuclideanDistance(Point3D point)
        {
            return Math.Sqrt(Math.Pow(X - point.X, 2) + Math.Pow(Y - point.Y, 2) + Math.Pow(Z - point.Z, 2));
        }

        public Point3D GetOppositeVector()
        {
            return new(-X, -Y, -Z);
        }

        public Point3D normalize()
        {

            double length = CalculateLength();
            if (length == 0)
            {
                return new(0, 0, 0);
            }
            return new(X / length, Y / length, Z / length);
        }

        public double CalculateScalarProduct(Point3D point)
        {
            return X * point.X + Y * point.Y + Z * point.Z;
        }

        public double CalculateLength()
        {
            return Math.Sqrt(X*X + Y*Y + Z*Z);
        }

        public double CalculateShortestDistanceToPlane(Plane2D plane)
        {
            Function planeEquation = plane.equation;
            double a = planeEquation.A;
            double b = planeEquation.B;
            double c = planeEquation.C;
            double d = planeEquation.D;
            return Math.Abs(a * X + b * Y + c * Z + d) / Math.Sqrt(a * a + b * b + c * c);

        }

        public Point3D CalculateCrossProduct(Point3D point)
        {
            return new(Y * point.Z - Z * point.Y,
                Z * point.X - X * point.Z,
                X * point.Y - Y * point.X);
        }

        public double CalculateDistanceFromOriginOfCoordinateSystem()
        {
            return Math.Sqrt((X * X) + (Y * Y) + (Z * Z));
        }

        public Point3D Multiply(double multiplicator)
        {
            return new(X * multiplicator, Y * multiplicator, Z * multiplicator);
        }

        public Point3D Subtract(Point3D point)
        {
            return new(X - point.X, Y - point.Y, Z - point.Z);
        }

        public Point3D Add(Point3D point)
        {
            return new(X + point.X, Y + point.Y, Z + point.Z);
        }

        public override bool Equals(object? obj)
        {
            return obj is Point3D d &&
                   X == d.X &&
                   Y == d.Y &&
                   Z == d.Z;
        }

        public override int GetHashCode()
        {
            return HashCode.Combine(X, Y, Z);
        }

        public override string? ToString()
        {
            return "Point3D(" + X + ", " + Y + ", " + Z + ")";
        }
    }
}
