using VoronoiDiagrams;

namespace CVD
{
    internal class VoronoiDiagram
    {

        public Dictionary<VoronoiPoint2D, VoronoiCell> CreateVoronoiDiagram(ISet<DelaunayTriangle> delaunayTriangulation)
        {
            Dictionary<Edge, List<DelaunayTriangle>> adjacencyTriangleMap = CreateAdjacencyTriangleMap(delaunayTriangulation);
            Dictionary<VoronoiPoint2D, VoronoiCell> voronoiCells = [];
  
            foreach (List<DelaunayTriangle> trianglesWithThisEdge in adjacencyTriangleMap.Values)
            {
                bool isEdgeSharedByTwoTriangles = trianglesWithThisEdge.Count == 2;
                if (isEdgeSharedByTwoTriangles)
                {
                    DelaunayTriangle triangle1 = trianglesWithThisEdge[0];
                    DelaunayTriangle triangle2 = trianglesWithThisEdge[1];

                    Edge edge1 = new(triangle1.circumCircleCenter, triangle2.circumCircleCenter);

                    AddVerticesIfAbsent(voronoiCells, triangle1);
                    AddVerticesIfAbsent(voronoiCells, triangle2);

                    VoronoiCell cell1 = voronoiCells[triangle1.point1];
                    VoronoiCell cell2 = voronoiCells[triangle1.point2];
                    VoronoiCell cell3 = voronoiCells[triangle1.point3];

                    cell1.AddEdge(edge1);
                    cell2.AddEdge(edge1);
                    cell3.AddEdge(edge1);
                }
            }

            return voronoiCells;
        }

        private static void AddVerticesIfAbsent(Dictionary<VoronoiPoint2D, VoronoiCell> voronoiCells, DelaunayTriangle triangle)
        {
            voronoiCells.TryAdd(triangle.point1, new(triangle.circumCircleCenter));
            voronoiCells.TryAdd(triangle.point2, new(triangle.circumCircleCenter));
            voronoiCells.TryAdd(triangle.point3, new(triangle.circumCircleCenter));
        }

        private Dictionary<Edge, List<DelaunayTriangle>> CreateAdjacencyTriangleMap(ISet<DelaunayTriangle> delaunayTriangulation)
        {
            Dictionary<Edge, List<DelaunayTriangle>> adjacencyTriangleMap = new();
            foreach (DelaunayTriangle triangle in delaunayTriangulation)
            {
                Edge edge1 = new Edge(triangle.point1, triangle.point2);
                PutToAdjacencyTrianglesMap(adjacencyTriangleMap, edge1, triangle);
                Edge edge2 = new Edge(triangle.point2, triangle.point3);
                PutToAdjacencyTrianglesMap(adjacencyTriangleMap, edge2, triangle);
                Edge edge3 = new Edge(triangle.point1, triangle.point3);
                PutToAdjacencyTrianglesMap(adjacencyTriangleMap, edge3 , triangle);

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
