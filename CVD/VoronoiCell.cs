using VoronoiDiagrams;

namespace CVD
{
    internal class VoronoiCell
    {
        private readonly List<Edge> edges = new();
        private readonly ISet<VoronoiPoint2D> vertices = new HashSet<VoronoiPoint2D>();

        private readonly VoronoiPoint2D center;


        public VoronoiCell(VoronoiPoint2D center)
        {
            this.center = center;
        }

        public void AddEdge(Edge edge)
        {
            edges.Add(edge);
            vertices.Add(edge.startingPoint);
            vertices.Add(edge.endingPoint);
        }

        public ISet<VoronoiPoint2D> getVertices() { return vertices; }

        public override bool Equals(object? obj)
        {
            return obj is VoronoiCell cell &&
                   EqualityComparer<List<Edge>>.Default.Equals(edges, cell.edges) &&
                   EqualityComparer<VoronoiPoint2D>.Default.Equals(center, cell.center);
        }

        public VoronoiPoint2D getCenter()
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
