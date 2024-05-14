using CVD.delaunay;
using CVD.shape;

namespace CVD.voronoi
{
    internal class VoronoiDiagram
    {

        public Dictionary<VoronoiPoint2D, VoronoiCell3D> CreateVoronoiDiagram(ISet<DelaunayTriangle> delaunayTriangulation)
        {
            Dictionary<Edge, List<DelaunayTriangle>> adjacencyTriangleMap = CreateAdjacencyTriangleMap(delaunayTriangulation);
            Dictionary<VoronoiPoint2D, VoronoiCell3D> voronoiCells = new();

            foreach (KeyValuePair<Edge, List<DelaunayTriangle>> keyPair in adjacencyTriangleMap)
            {
                Edge currentEdge = keyPair.Key;
                List<DelaunayTriangle> trianglesWithThisEdge = keyPair.Value;
                bool isEdgeSharedByTwoTriangles = trianglesWithThisEdge.Count == 2;
                if (isEdgeSharedByTwoTriangles)
                {
                    DelaunayTriangle triangle1 = trianglesWithThisEdge[0];
                    DelaunayTriangle triangle2 = trianglesWithThisEdge[1];


                    Edge voronoiCellEdge = new(triangle1.circumCircleCenter, triangle2.circumCircleCenter);

                    AddVerticesIfAbsent(voronoiCells, currentEdge);

                    VoronoiCell3D cell1 = voronoiCells[currentEdge.startingPoint];
                    VoronoiCell3D cell2 = voronoiCells[currentEdge.endingPoint];

                    Edge3D voronoiCellEdge3D = CreateEdge3D(voronoiCellEdge);

                    cell1.AddEdge(voronoiCellEdge3D);
                    cell2.AddEdge(voronoiCellEdge3D);
                }
            }

            return voronoiCells;
        }

        private Edge3D CreateEdge3D(Edge edge)
        {
            VoronoiPoint2D startingPoint = edge.startingPoint;
            VoronoiPoint2D endingPoint = edge.endingPoint;
            return new(new(startingPoint.X, startingPoint.Y, 0), new(endingPoint.X, endingPoint.Y, 0));
        }

        private static void AddVerticesIfAbsent(Dictionary<VoronoiPoint2D, VoronoiCell3D> voronoiCells, Edge edge)
        {
            Point3D startingPoint = new(edge.startingPoint.X, edge.startingPoint.Y, 0);
            Point3D endingPoint = new(edge.endingPoint.X, edge.endingPoint.Y, 0);

            voronoiCells.TryAdd(edge.startingPoint, new(startingPoint));
            voronoiCells.TryAdd(edge.endingPoint, new(endingPoint));
        }

        private Dictionary<Edge, List<DelaunayTriangle>> CreateAdjacencyTriangleMap(ISet<DelaunayTriangle> delaunayTriangulation)
        {
            Dictionary<Edge, List<DelaunayTriangle>> adjacencyTriangleMap = new();
            foreach (DelaunayTriangle triangle in delaunayTriangulation)
            {
                Edge edge1 = new(triangle.point1, triangle.point2);
                PutToAdjacencyTrianglesMap(adjacencyTriangleMap, edge1, triangle);
                Edge edge2 = new(triangle.point2, triangle.point3);
                PutToAdjacencyTrianglesMap(adjacencyTriangleMap, edge2, triangle);
                Edge edge3 = new(triangle.point1, triangle.point3);
                PutToAdjacencyTrianglesMap(adjacencyTriangleMap, edge3, triangle);
            }

            return adjacencyTriangleMap;
        }

        private static void PutToAdjacencyTrianglesMap(Dictionary<Edge, List<DelaunayTriangle>> map, Edge edge,
                                                   DelaunayTriangle delaunayTriangle)
        {
            if (map.ContainsKey(edge))
            {
                map[edge].Add(delaunayTriangle);
                return;
            }
            List<DelaunayTriangle> value = new();
            value.Add(delaunayTriangle);
            map[edge] = value;
        }
    }
}
