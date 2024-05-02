using VoronoiDiagrams;

namespace CVD
{
    internal class VoronoiDiagram
    {

        public Dictionary<Point3D, VoronoiCell> CreateVoronoiDiagram(ISet<DelaunayTriangle> delaunayTriangulation)
        {
            Dictionary<Edge, List<DelaunayTriangle>> adjacencyTriangleMap = CreateAdjacencyTriangleMap(delaunayTriangulation);
            Dictionary<Point3D, VoronoiCell> voronoiCells = new();
  
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

                    VoronoiCell cell1 = voronoiCells[currentEdge.startingPoint];
                    VoronoiCell cell2 = voronoiCells[currentEdge.endingPoint];
                    
                    cell1.AddEdge(voronoiCellEdge);
                    cell2.AddEdge(voronoiCellEdge);             
                }
            }

            return voronoiCells;
        }

        private static void AddVerticesIfAbsent(Dictionary<Point3D, VoronoiCell> voronoiCells, Edge edge)
        {
            voronoiCells.TryAdd(edge.startingPoint, new(edge.startingPoint));
            voronoiCells.TryAdd(edge.endingPoint, new(edge.endingPoint));
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
