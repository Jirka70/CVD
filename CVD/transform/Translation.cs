using CVD.shape;

namespace CVD.transforms
{
    internal class Translation : Transform, RevertableTransform
    {

        private Point3D? _translation = null;


        public Point3D Transform(Point3D point)
        {
            if (_translation == null) // debilni, trochu predelat...
            {
                return new(point.X, point.Y, 0);
            }

            return new(point.X + _translation.X, point.Y + _translation.Y, point.Z + _translation.Z);
        }

        public void SetTranslation(Point3D translation)
        {
            _translation = translation;
        }

        public Point3D RevertTransform(Point3D point)
        {
            if (_translation == null)
            {
                return point;
            }

            return new(point.X - _translation.X, point.Y - _translation.Y, point.Z - _translation.Z);
        }
    }
}
