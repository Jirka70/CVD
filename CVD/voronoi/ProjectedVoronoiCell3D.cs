using CVD.shape;

namespace CVD.voronoi
{
    internal class ProjectedVoronoiCell3D
    {
        public readonly Point3D center;
        public readonly ISet<Edge3D> edgesInRectangle = new HashSet<Edge3D>();
        public readonly ISet<Edge3D> edgesOutRectangle = new HashSet<Edge3D>();

        public ProjectedVoronoiCell3D(Point3D center)
        {
            this.center = center; 
        }

        public void AddInRectangleEdge(Edge3D edge)
        {
            edgesInRectangle.Add(edge);
        }

        public void AddOutRectangleEdge(Edge3D edge)
        {
            edgesOutRectangle.Add(edge);
        }
    }
}
