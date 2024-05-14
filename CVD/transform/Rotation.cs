using CVD.shape;

namespace CVD.transforms
{
    internal class Rotation : Transform, RevertableTransform
    {
        private readonly Point3D _rotationAxis;
        private readonly Point3D[] _rotationMatrix;
        private readonly Point3D[] _negativeAngleRotationMatrix;

        private float _rotationAngle = 0f;
        private float _negativeRotationAngle;

        public Rotation(Point3D rotationAxis, float rotationAngleRads)
        {
            _rotationAxis = rotationAxis.normalize();
            _rotationAngle = rotationAngleRads;
            _negativeRotationAngle = -_rotationAngle;

            double x = _rotationAxis.X;
            double y = _rotationAxis.Y;
            double z = _rotationAxis.Z;

            double cosAngle = Math.Cos(_rotationAngle);
            double sinAngle = Math.Sin(_rotationAngle);

            _rotationMatrix = new Point3D[] {
                new(cosAngle + x * x * (1 - cosAngle), x * y * (1 - cosAngle) - z * sinAngle, x * z * (1 - cosAngle) + y * sinAngle),
                new(y * x * (1 - cosAngle) + z * sinAngle, cosAngle + y * y * (1 - cosAngle), y * z * (1 - cosAngle) - x * sinAngle),
                new(z * x * (1 - cosAngle) - y * sinAngle, z * y * (1 - cosAngle) + x * sinAngle, cosAngle + z * z * (1 - cosAngle))};

            double negativeCosAngle = Math.Cos(_negativeRotationAngle);
            double negativeSinAngle = Math.Sin(_negativeRotationAngle);

            _negativeAngleRotationMatrix = new Point3D[]
            {
                new(negativeCosAngle + x * x * (1 - negativeCosAngle), x * y * (1 - negativeCosAngle) - z * negativeSinAngle, x * z * (1 - negativeCosAngle) + y * negativeSinAngle),
                new(y * x * (1 - negativeCosAngle) + z * negativeSinAngle, negativeCosAngle + y * y * (1 - negativeCosAngle), y * z * (1 - negativeCosAngle) - x * negativeSinAngle),
                new(z * x * (1 - negativeCosAngle) - y * negativeSinAngle, z * y * (1 - negativeCosAngle) + x * negativeSinAngle, negativeCosAngle + z * z * (1 - negativeCosAngle))};
        }

        public Point3D RevertTransform(Point3D point)
        {
            double first = _negativeAngleRotationMatrix[0].CalculateScalarProduct(point);
            double second = _negativeAngleRotationMatrix[1].CalculateScalarProduct(point);
            double third = _negativeAngleRotationMatrix[2].CalculateScalarProduct(point);

            return new(first, second, third);
            throw new NotImplementedException();
        }

        public Point3D Transform(Point3D point)
        {
            double first = _rotationMatrix[0].CalculateScalarProduct(point);
            double second = _rotationMatrix[1].CalculateScalarProduct(point);
            double third = _rotationMatrix[2].CalculateScalarProduct(point);

            return new(first, second, third);
        }
    }
}
