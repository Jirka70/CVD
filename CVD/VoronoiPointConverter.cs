namespace VoronoiDiagrams
{
    internal class VoronoiPointConverter
    {
        public static List<VoronoiPoint2D> Convert3DPointsTo2DPoints(List<Point3D> points3D)
        {
            List<VoronoiPoint2D> convertedPoints = new();
            foreach (Point3D point3D in points3D)
            {
                VoronoiPoint2D convertedPoint = Convert3DPointTo2DPoint(point3D);
                convertedPoints.Add(convertedPoint);
            }

            return convertedPoints;
        }

        /**
         * Z coord will be completely neglected
         */
        private static VoronoiPoint2D Convert3DPointTo2DPoint(Point3D point3D)
        {
            return new(point3D.x, point3D.y);
        }
    }
}
