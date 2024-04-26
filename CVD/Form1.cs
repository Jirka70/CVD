using System.Collections.Generic;
using System.Drawing;
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
            List<VoronoiPoint2D> projectedPointCloud = VoronoiPointConverter.Convert3DPointsTo2DPoints(randomPointCloud);
       
            DelaunayTriangulation delaunayTriangulation = new();
            ISet<DelaunayTriangle> triangulation = delaunayTriangulation.CreateTriangulation(projectedPointCloud);
            VoronoiDiagram voronoi = new();
            Dictionary<VoronoiPoint2D, VoronoiCell> voronoiDiagram = voronoi.CreateVoronoiDiagram(triangulation);

            voronoiDiagram = CenterVoronoiDiagram(voronoiDiagram, 42);
            DrawVoronoiCells(voronoiDiagram);


            /*List <Point3D> randomPointCloud = Point3DGenerator.GenerateRandomPoints(POINT_CLOUD_SIZE, WIDTH, HEIGHT);

            List<VoronoiPoint2D> points = new();
            for (int i = 0; i < HEIGHT; i += 50)
            {
                for (int j = 0; j < WIDTH; j += 50)
                {
                    points.Add(new(j, i));
                }
            }
            Func<double, double> slopeFunction = x => 2*x;
            List<VoronoiPoint2D> transformedPoints = TransformPointsByLinearFunction(points, slopeFunction);
            DelaunayTriangulation delaunayTriangulation = new();
            delaunayTriangulation.CreateTriangulation(transformedPoints);
            Draw3DPoints(points);  */ 
        }

        private void Draw3DPoints(List<VoronoiPoint2D> points)
        {
            Brush brush = new SolidBrush(Color.Black);
            foreach (VoronoiPoint2D point in points)
            {
                Draw3DPoints(brush, point);
            }
        }

        private void Draw3DPoints(Brush brush, VoronoiPoint2D point)
        {
            int pointWidth = 4;
            int pointHeight = 4;
            if (point.X > 0 && point.Y > 0)
            {
                graphicsContext.FillEllipse(brush, (int)point.X - pointWidth / 2, (int)point.Y - pointHeight / 2, pointWidth, pointHeight);
            }
        }

        private List<VoronoiPoint2D> TransformPointsByLinearFunction(List<VoronoiPoint2D> points, Func<double, double> slopeFunction)
        {
            List<VoronoiPoint2D> transformedPoints = new();
            
            /*for (int i = 0; i < points.Count; i++)
            {
                VoronoiPoint2D currentPoint = points[i];

                double x = slopeFunction(currentPoint.X);
                double derivation = x / currentPoint.X;
                x = currentPoint.X / derivation;
                double y = currentPoint.Y;
                double z = slopeFunction(currentPoint.X);

                VoronoiPoint2D transformedPoint = new(x, y, z);
                transformedPoints.Add(transformedPoint);
            }*/

            return transformedPoints;
        }

        private Dictionary<VoronoiPoint2D, VoronoiCell> CalculateVoronoiDiagram(List<VoronoiPoint2D> points)
        {
            DelaunayTriangulation delaunayTriangulation = new();
            ISet<DelaunayTriangle> triangulation = delaunayTriangulation.CreateTriangulation(points);
            
            VoronoiDiagram voronoi = new();
            return voronoi.CreateVoronoiDiagram(triangulation);
        }

        private Dictionary<VoronoiPoint2D, VoronoiCell> CenterVoronoiDiagram(Dictionary<VoronoiPoint2D, VoronoiCell> voronoiDiagram, int numberOfIterations)
        {
            List<VoronoiPoint2D> newPoints = new();
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
            DrawPoints(newPoints);
            
           
            return voronoiDiagram;
        }

        private VoronoiPoint2D CalculateCenterOfGravity(VoronoiCell voronoiCell)
        {
            ISet<VoronoiPoint2D> vertices = voronoiCell.getVertices();
            int verticesSize = vertices.Count;
            double findingX = 0;
            double findingY = 0;
            foreach (VoronoiPoint2D vertex in vertices)
            {
                if (Math.Abs(vertex.X) > WIDTH || Math.Abs(vertex.Y) > HEIGHT)
                {
                    return voronoiCell.getCenter();
                }
                findingX += vertex.X;
                findingY += vertex.Y;
            }

            return new VoronoiPoint2D(findingX / verticesSize, findingY / verticesSize);
        }

        private void DrawVoronoiCells(Dictionary<VoronoiPoint2D, VoronoiCell> voronoiCells)
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
