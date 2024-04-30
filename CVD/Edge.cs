using VoronoiDiagrams;

namespace CVD
{
    internal class Edge
    {
        public readonly Point3D startingPoint;
        public readonly Point3D endingPoint;

        public Edge(Point3D startingPoint, Point3D endingPoint)
        {
            this.startingPoint = startingPoint;
            this.endingPoint = endingPoint;
        }

        public bool HasVertex(Point3D vertex)
        {
            return startingPoint.Equals(vertex) || endingPoint.Equals(vertex);
        }

        public override bool Equals(object? obj)
        {
            if (obj == null || GetType() != obj.GetType())
            {
                return false;
            }

            Edge other = (Edge)obj;

            return startingPoint.Equals(other.startingPoint) && endingPoint.Equals(other.endingPoint);
        }

        public override int GetHashCode()
        {
            return HashCode.Combine(startingPoint, endingPoint);
        }
    }
}
