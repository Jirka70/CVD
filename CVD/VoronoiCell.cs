using VoronoiDiagrams;

namespace CVD
{
    internal class VoronoiCell
    {
        private readonly List<Edge> edges = new();
        private readonly ISet<Point3D> vertices = new HashSet<Point3D>();

        private readonly Point3D center;


        public VoronoiCell(Point3D center)
        {
            this.center = center;
        }

        public void AddEdge(Edge edge)
        {
            edges.Add(edge);
            vertices.Add(edge.startingPoint);
            vertices.Add(edge.endingPoint);
        }

        public ISet<Point3D> getVertices() { return vertices; }

        public override bool Equals(object? obj)
        {
            return obj is VoronoiCell cell &&
                   EqualityComparer<List<Edge>>.Default.Equals(edges, cell.edges) &&
                   EqualityComparer<Point3D>.Default.Equals(center, cell.center);
        }

        public Point3D getCenter()
        {
            return center;
        }

        public override int GetHashCode()
        {
            return HashCode.Combine(vertices, center);
        }

        public List<Edge> GetEdges() 
        { 
            return edges;
        }

        
    }
}
