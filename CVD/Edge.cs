using VoronoiDiagrams;

namespace CVD
{
    internal class Edge
    {
        public readonly VoronoiPoint2D startingPoint;
        public readonly VoronoiPoint2D endingPoint;

        public Edge(VoronoiPoint2D startingPoint, VoronoiPoint2D endingPoint)
        {
            this.startingPoint = startingPoint;
            this.endingPoint = endingPoint;
        }

        public bool HasVertex(VoronoiPoint2D vertex)
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
