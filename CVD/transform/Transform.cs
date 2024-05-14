using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using CVD.shape;

namespace CVD.transforms
{
    internal interface Transform
    {
        Point3D Transform(Point3D point);
    }
}
