using VoronoiDiagrams;

namespace CVD
{
    public partial class VoronoiForm : Form
    {
        private static readonly int POINT_CLOUD_SIZE = 500;

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
            List<Point3D> randomPointCloud = Point3DGenerator.GenerateRandomPoints(POINT_CLOUD_SIZE, WIDTH, HEIGHT);
            /*List<Point3D> randomPointCloud = new();
            for (int i = 0; i < HEIGHT; i += 50)
            {
                for (int j = 0; j < WIDTH; j += 50)
                {
                    randomPointCloud.Add(new(j, i, 12));
                }
            }*/
            List<Point3D> projectedPointCloud = CalculatePointsOnSlope(0, randomPointCloud);
            DelaunayTriangulation delaunayTriangulation = new();
            ISet<DelaunayTriangle> triangulation = delaunayTriangulation.CreateTriangulation(projectedPointCloud);
            //DrawTriangles(triangulation);
            VoronoiDiagram voronoi = new();
            Dictionary<Point3D, VoronoiCell> voronoiDiagram = voronoi.CreateVoronoiDiagram(triangulation);
            voronoiDiagram = CenterVoronoiDiagram(voronoiDiagram, 50);
            DrawVoronoiCells(voronoiDiagram);
            //Draw3DPoints(projectedPointCloud);
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
