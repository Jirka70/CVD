using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using VoronoiDiagrams;

namespace CVD
{
    internal class Plane2D
    {
        public readonly Function equation;
        public readonly Point3D normalVector;
        public Plane2D(Point3D normalVector, Point3D initialPoint)
        {
            this.normalVector = normalVector;
            double a = normalVector.X;
            double b = normalVector.Y;
            double c = normalVector.Z;
            double d = -(initialPoint.X * a + initialPoint.Y * b + initialPoint.Z * c);
            equation = new Function(a, b, c, d);
        }

        public bool ContainsPoint(Point3D point)
        {
            double epsilon = 0.0001;
            double pointValue = equation.CalculateValueOfPoint(point);
            return Math.Abs(pointValue) <= epsilon;        
        }

        /**
         * @return smallest angle between two planes in RADIANS          
         */
        public double CalculateAngle(Plane2D plane)
        {
            double scalarProduct = normalVector.CalculateScalarProduct(plane.normalVector);
            double normalVectorLength = normalVector.CalculateLength();
            double otherNormalVectorLength = plane.normalVector.CalculateLength();
            return Math.Acos(scalarProduct / (normalVectorLength * otherNormalVectorLength));
        }
    }
}
