﻿using VoronoiDiagrams;

namespace CVD
{
    internal class DelaunayTriangulation
    {
        private static readonly int BOUND = 1_000_000;
        private readonly DelaunayTriangle superDelaunayTriangle;

        public DelaunayTriangulation() {
            superDelaunayTriangle = new(new Point3D(BOUND / 2, -BOUND, 0), new Point3D(-BOUND, BOUND, 0),
                new Point3D(BOUND, BOUND, 0));   
        }

        public ISet<DelaunayTriangle> CreateTriangulation(List<Point3D> points)
        {
            ISet<DelaunayTriangle> delaunayTriangles = new HashSet<DelaunayTriangle>();
            delaunayTriangles.Add(superDelaunayTriangle);
            foreach (Point3D point in points)
            {
                ISet<DelaunayTriangle> badTriangles = FindBadTriangles(point, delaunayTriangles);
                ISet<Edge> polygon = CreatePolygonOfBadTriangles(badTriangles);
                RemoveBadTrianglesFromTriangulation(delaunayTriangles, badTriangles);
                Triangulate(point, delaunayTriangles, polygon);
            }

            return delaunayTriangles;
        }

        private void RemoveTrianglesContainingSupertriangleVertex(ISet<DelaunayTriangle> triangles)
        {
            Point3D superTriangleVertex1 = superDelaunayTriangle.point1;
            Point3D superTriangleVertex2 = superDelaunayTriangle.point2;
            Point3D superTriangleVertex3 = superDelaunayTriangle.point3;

            foreach (DelaunayTriangle triangle in triangles) {
                if (triangle.HasVertex(superTriangleVertex1)
                    || triangle.HasVertex(superTriangleVertex2)
                    || triangle.HasVertex(superTriangleVertex3))
                {
                    triangles.Remove(triangle);
                }
            }
        }

        private static void Triangulate(Point3D point, ISet<DelaunayTriangle> triangulation, ISet<Edge> polygon)
        {
            foreach(Edge edge in polygon)
            {
                triangulation.Add(new DelaunayTriangle(edge.startingPoint, edge.endingPoint, new(point.X, point.Y, point.Z)));
            }
        }

        private static void RemoveBadTrianglesFromTriangulation(ISet<DelaunayTriangle> delaunayTriangles,
                                                     ISet<DelaunayTriangle> badDelaunayTriangles)
        {
            delaunayTriangles.ExceptWith(badDelaunayTriangles);
        }

        private static bool IsEdgeSharedWithOtherBadTriangle(ISet<DelaunayTriangle> badTriangles,
                                                     DelaunayTriangle triangleWithEdge, Edge edge)
        {
            foreach(DelaunayTriangle badTriangle in badTriangles)
            {
                if (badTriangle.Equals(triangleWithEdge))
                {
                    continue;
                }
                if (badTriangle.HasEdge(edge))
                {
                    return true;
                }
            }
            return false;
        }



        private static ISet<Edge> CreatePolygonOfBadTriangles(ISet<DelaunayTriangle> badTriangles)
        {
            ISet<Edge> polygon = new HashSet<Edge>();
            foreach (DelaunayTriangle badTriangle in badTriangles)
            {

                Edge edge1 = new(badTriangle.point1, badTriangle.point2);
                Edge edge2 = new(badTriangle.point2, badTriangle.point3);
                Edge edge3 = new(badTriangle.point1, badTriangle.point3);
                if (!IsEdgeSharedWithOtherBadTriangle(badTriangles, badTriangle, edge1))
                {
                    polygon.Add(edge1);
                }
                if (!IsEdgeSharedWithOtherBadTriangle(badTriangles, badTriangle, edge2))
                {
                    polygon.Add(edge2);
                }
                if (!IsEdgeSharedWithOtherBadTriangle(badTriangles, badTriangle, edge3))
                {
                    polygon.Add(edge3); 
                }
            }

            return polygon;
        }

        private static ISet<DelaunayTriangle> FindBadTriangles(Point3D point, ISet<DelaunayTriangle> triangulation)
        {
            ISet<DelaunayTriangle> badTriangles = new HashSet<DelaunayTriangle>();
            foreach (DelaunayTriangle triangle in triangulation)
            {
                if (triangle.IsInCircumcircle(point.X, point.Y, point.Z))
                {
                    badTriangles.Add(triangle);
                }
            }

            return badTriangles;
        }
    }
}
