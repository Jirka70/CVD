using CVD.shape;
using CVD.transforms;

namespace CVD.transform
{
    internal class PointTransformation : Transform, RevertableTransform
    {
        private readonly RevertableTransform[] _transforms;

        public PointTransformation(params RevertableTransform[] transforms)
        {
            _transforms = transforms;
        }

        public Point3D RevertTransform(Point3D point)
        {
            foreach (RevertableTransform transform in _transforms.Reverse())
            {
                point = transform.RevertTransform(point);
            }

            return point;
        }

        public Point3D Transform(Point3D point)
        {
            foreach (Transform transform in _transforms)
            {
                point = transform.Transform(point);
            }

            return point;
        }
    }
}
