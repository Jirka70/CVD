using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CVD.shape
{
    internal class Edge3D
    {
        public readonly Point3D startingPoint;
        public readonly Point3D endingPoint;
        public Edge3D(Point3D startingPoint, Point3D endingPoint)
        {
            this.startingPoint = startingPoint;
            this.endingPoint = endingPoint;
        }
    }
}
