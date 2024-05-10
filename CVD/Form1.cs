using System.Numerics;
using VoronoiDiagrams;

namespace CVD
{
    public partial class VoronoiForm : Form
    {
        private static readonly int POINT_CLOUD_SIZE = 42;

        private readonly Graphics graphicsContext;
        private readonly PictureBox canvasBox = new PictureBox();


        public VoronoiForm()
        {
            InitializeComponent();
            canvasBox.Image = new Bitmap(WIDTH, HEIGHT);
            canvasBox.Size = new Size(WIDTH, HEIGHT);
            canvasBox.Location = new Point(0,0);
            Controls.Add(canvasBox);
            Bitmap image = (Bitmap) canvasBox.Image;
            graphicsContext = Graphics.FromImage(image);

            VoronoiProcess();
        }



        private void VoronoiProcess()
        {
            List<Point3D> randomPointCloud = new();//Point3DGenerator.GenerateRandomPoints(POINT_CLOUD_SIZE, WIDTH, HEIGHT);
            Point3D pointToProject = new(3,3,3);
            randomPointCloud.Add(pointToProject);
            Point3D initialPoint = new(1,2,3);
            Point3D normalVector = new(0,0,1);
            Point3D oppositeVector = normalVector.GetOppositeVector();
            Plane2D plane2D = new(normalVector, initialPoint);
            MessageBox.Show(plane2D.equation.ToString());

            /*double distanceToPlane = point.CalculateShortestDistanceToPlane(plane2D);
            Point3D oppositeVectorNormalized = oppositeVector.normalize();
            Point3D projectedPoint = oppositeVectorNormalized.Multiply(distanceToPlane).Add(point);
            if (!plane2D.ContainsPoint(projectedPoint))
            {
                projectedPoint = normalVector.normalize().Multiply(distanceToPlane).Add(point);
            }
            
            Point3D pedestalPlaneNormalizedVector = new(0, 0, 1);

            Point3D axisVector = normalVector.CalculateCrossProduct(pedestalPlaneNormalizedVector);
            Point3D intersectPoint = CalculateIntersectPointWithPedestalPlane(plane2D);
            Point3D translatedProjectedPoint = projectedPoint.Subtract(intersectPoint);
            Plane2D pedestalPlane = new(new(0, 0, 1), new(0, 0, 0));
            double angle = plane2D.CalculateAngle(pedestalPlane);
            Point3D rotatedPoint = RotatePointAroundArbitraryAxis(axisVector, translatedProjectedPoint, angle);
            MessageBox.Show("" + rotatedPoint.ToString());*/
            List<Point3D> projectedPointloud = ProjectPointCloudOnGround(plane2D, randomPointCloud);



            /*List<Point3D> randomPointCloud = new();
            for (int i = 0; i < HEIGHT; i += 50)
            {
                for (int j = 0; j < WIDTH; j += 50)
                {
                    randomPointCloud.Add(new(j, i, 12));
                }
            }*/
            /*List<Point3D> projectedPointCloud = CalculatePointsOnSlope(0, randomPointCloud);
            DelaunayTriangulation delaunayTriangulation = new();
            ISet<DelaunayTriangle> triangulation = delaunayTriangulation.CreateTriangulation(projectedPointCloud);
            //DrawTriangles(triangulation);
            VoronoiDiagram voronoi = new();
            Dictionary<Point3D, VoronoiCell> voronoiDiagram = voronoi.CreateVoronoiDiagram(triangulation);
            voronoiDiagram = CenterVoronoiDiagram(voronoiDiagram, 50);
            DrawVoronoiCells(voronoiDiagram);
            //Draw3DPoints(projectedPointCloud);*/
        }

        private List<Point3D> ProjectPointCloudOnGround(Plane2D plane, List<Point3D> pointcloud)
        {
            List<Point3D> projectedPointcloud = new();
            foreach (Point3D point in pointcloud)
            {
                Point3D projectedPoint = ProjectPointOnGround(plane, point);
                projectedPointcloud.Add(projectedPoint);
                MessageBox.Show("pp: " + projectedPoint.ToString());
            }

            return projectedPointcloud;
        }



        private Point3D ProjectPointOnGround(Plane2D plane, Point3D point)
        {
            Point3D normalVector = plane.normalVector;
            Point3D oppositeVector = normalVector.GetOppositeVector();
            double distanceToPlane = point.CalculateShortestDistanceToPlane(plane);
            Point3D oppositeVectorNormalized = oppositeVector.normalize();
            Point3D projectedPoint = oppositeVectorNormalized.Multiply(distanceToPlane).Add(point);

            if (!plane.ContainsPoint(projectedPoint))
            {
                projectedPoint = normalVector.normalize().Multiply(distanceToPlane).Add(point);
            }

            Point3D pedestalPlaneNormalizedVector = new(0, 0, 1);
            Point3D axisVector = normalVector.CalculateCrossProduct(pedestalPlaneNormalizedVector);
            Point3D intersectPoint = CalculateIntersectPointWithPedestalPlane(plane);
            Point3D translatedProjectedPoint = projectedPoint.Subtract(intersectPoint);

            Plane2D pedestalPlane = new(new(0, 0, 1), new(0, 0, 0));
            double angle = plane.CalculateAngle(pedestalPlane);
            return RotatePointAroundArbitraryAxis(axisVector, translatedProjectedPoint, angle);
        } 

        private Point3D CalculateIntersectPointWithPedestalPlane(Plane2D plane)
        {
            double a = plane.equation.A;
            double b = plane.equation.B;
            double d = plane.equation.D;
            if (a == 0 && b == 0)
            {
                return new(0,0,0);
            }
            if (b == 0)
            {
                return new(-d / a, 0, 0);
            }
            double x = 2;
            double y = (-d - a*x) / b;
            return new(x, y, 0);
        }

        /**
         * Angle has to be defined in RADIANS
         */
        private Point3D RotatePointAroundArbitraryAxis(Point3D axis, Point3D point, double angle)
        {
            Point3D axisNormalized = axis.normalize();
            double x = axisNormalized.X;
            double y = axisNormalized.Y;
            double z = axisNormalized.Z;
            double cosAngle = Math.Cos(angle);
            double sinAngle = Math.Sin(angle);
            Point3D[] rotationMatrix = { new(cosAngle + x*x*(1 - cosAngle), x*y*(1-cosAngle) - z*sinAngle, x*z*(1-cosAngle) + y*sinAngle),
                new(y*x*(1-cosAngle) + z*(sinAngle), cosAngle+y*y*(1-cosAngle), y*z*(1-cosAngle) - x*sinAngle),
                new(z*x*(1-cosAngle)-y*sinAngle, z*y*(1-cosAngle)+x*sinAngle, cosAngle + z*z*(1-cosAngle))};
            return Rotate(rotationMatrix, point);
        }

        private Point3D Rotate(Point3D[] rotationMatrix, Point3D vectorToRotate)
        {
            double first = rotationMatrix[0].CalculateScalarProduct(vectorToRotate);
            double second = rotationMatrix[1].CalculateScalarProduct(vectorToRotate);
            double third = rotationMatrix[2].CalculateScalarProduct(vectorToRotate);
            

            return new(first, second, third);
        }



        private void Draw3DPoints(List<Point3D> points)
        {
            Brush brush = new SolidBrush(Color.Black);
            foreach (Point3D point in points)
            {
                Draw3DPoints(brush, point);
            }
        }

        private void Draw3DPoints(Brush brush, Point3D point)
        {
            int pointWidth = 4;
            int pointHeight = 4;
            if (point.X > 0 && point.Y > 0)
            {
                graphicsContext.FillEllipse(brush, (int)point.X - pointWidth / 2, (int)point.Y - pointHeight / 2, pointWidth, pointHeight);
            }
        }

        private List<Point3D> CalculatePointsOnSlope(double slopeDegreesAngle, List<Point3D> points)
        {
            double angle = slopeDegreesAngle * (Math.PI / 180);

            List<Point3D> newPoints = new();
            foreach (Point3D point in points)
            {
                double x = point.X;
                double y = point.Y;

                double newX = x * Math.Cos(angle);
                double newY = y;
                double newZ = x * Math.Tan(angle);
                Point3D newPoint = new(newX, newY, newZ);
                newPoints.Add(newPoint);
            }

            return newPoints;

        }



        private Dictionary<Point3D, VoronoiCell> CalculateVoronoiDiagram(List<Point3D> points)
        {
            DelaunayTriangulation delaunayTriangulation = new();
            ISet<DelaunayTriangle> triangulation = delaunayTriangulation.CreateTriangulation(points);
            
            VoronoiDiagram voronoi = new();
            return voronoi.CreateVoronoiDiagram(triangulation);
        }

        private Dictionary<Point3D, VoronoiCell> CenterVoronoiDiagram(Dictionary<Point3D, VoronoiCell> voronoiDiagram, int numberOfIterations)
        {
            List<Point3D> newPoints = new();
            for (int i = 0; i < numberOfIterations; i++)
            {
                newPoints = new();
                foreach (VoronoiCell cell in voronoiDiagram.Values)
                {
                    newPoints.Add(CalculateCenterOfGravity(cell));
                }
                voronoiDiagram = CalculateVoronoiDiagram(newPoints);
            }
            DrawVoronoiCells(voronoiDiagram);
            Draw3DPoints(newPoints);
            
           
            return voronoiDiagram;
        }

        private Point3D CalculateCenterOfGravity(VoronoiCell voronoiCell)
        {
            ISet<Point3D> vertices = voronoiCell.getVertices();
            int verticesSize = vertices.Count;
            double findingX = 0;
            double findingY = 0;
            double findingZ = 0;
            foreach (Point3D vertex in vertices)
            {
                if (Math.Abs(vertex.X) > WIDTH || Math.Abs(vertex.Y) > HEIGHT)
                {
                    return voronoiCell.getCenter();
                }
                findingX += vertex.X;
                findingY += vertex.Y;
                findingZ += vertex.Z;
            }

            return new(findingX / verticesSize, findingY / verticesSize, findingZ);
        }

        private void DrawVoronoiCells(Dictionary<Point3D, VoronoiCell> voronoiCells)
        {
            Pen pen = new(Color.Blue);
            foreach (VoronoiCell cell in voronoiCells.Values)
            {             
                DrawVoronoiCell(pen, cell);
            }
        }

        private void DrawVoronoiCell(Pen pen, VoronoiCell cell)
        {
            foreach (Edge vertex in cell.GetEdges())
            {
              
                graphicsContext.DrawLine(pen,
                    new((int) vertex.startingPoint.X, (int) vertex.startingPoint.Y),
                    new((int) vertex.endingPoint.X, (int) vertex.endingPoint.Y));
            }
        }

        private void DrawPoints(List<VoronoiPoint2D> points)
        {
            Brush brush = new SolidBrush(Color.Black);
            foreach(VoronoiPoint2D point in points) {
                DrawPoint(brush, point);
            }
        }

        private void DrawTriangles(ICollection<DelaunayTriangle> triangles)
        {
            Pen pen = new Pen(Color.Red);
            foreach(DelaunayTriangle triangle in triangles)
            {
                DrawTriangle(pen, triangle);
            }
        }

        private void DrawPoint(Brush brush, VoronoiPoint2D point)
        {

            int pointWidth = 4;
            int pointHeight = 4;
            if (point.X > 0 && point.Y > 0)
            {
                graphicsContext.FillEllipse(brush, (int)point.X - pointWidth / 2, (int)point.Y - pointHeight / 2, pointWidth, pointHeight);
            }
        }

        private void DrawTriangle(Pen pen, DelaunayTriangle triangle)
        {
            Point[] points = { new((int) triangle.point1.X, (int) triangle.point1.Y),
                new((int)triangle.point2.X, (int) triangle.point2.Y),
                new((int) triangle.point3.X, (int) triangle.point3.Y) };
            graphicsContext.DrawPolygon(pen, points);
        }
    }
}
