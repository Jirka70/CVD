using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using VoronoiDiagrams;

namespace CVD
{
    internal class Function
    {
        public readonly double A;
        public readonly double B;
        public readonly double C;
        public readonly double D;


        public Function(double a, double b, double c, double d)
        {
            A = a;
            B = b;
            C = c;
            D = d;
        }

        public double CalculateValueOfPoint(Point3D point)
        {
            return A*point.X + B*point.Y + C*point.Z + D;
        }

        public override string? ToString()
        {
            return A + "x + (" + B + ")y + (" + C + ")z + (" + D + ") = 0s";
        }
    }
}
