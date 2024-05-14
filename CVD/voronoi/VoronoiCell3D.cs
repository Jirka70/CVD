using CVD.shape;

namespace CVD.voronoi
{
    internal class VoronoiCell3D
    {
        private readonly List<Edge3D> edges = new();
        private readonly ISet<Point3D> vertices = new HashSet<Point3D>();

        private readonly Point3D center;


        public VoronoiCell3D(Point3D center)
        {
            this.center = center;
        }

        public void AddEdge(Edge3D edge)
        {
            edges.Add(edge);
            Point3D startingPoint = new(edge.startingPoint.X, edge.startingPoint.Y, 0);
            Point3D endingPoint = new(edge.endingPoint.X, edge.endingPoint.Y, 0);
            vertices.Add(startingPoint);
            vertices.Add(endingPoint);
        }

        public ISet<Point3D> getVertices() { return vertices; }

        public override bool Equals(object? obj)
        {
            return obj is VoronoiCell3D cell &&
                   EqualityComparer<List<Edge3D>>.Default.Equals(edges, cell.edges) &&
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

        public List<Edge3D> GetEdges()
        {
            return edges;
        }
    }
}
