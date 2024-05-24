using System.Diagnostics;
using System.Drawing;
using CVD.delaunay;
using CVD.shape;
using CVD.transform;
using CVD.transforms;
using CVD.util;
using CVD.voronoi;
using static System.Windows.Forms.DataFormats;

namespace CVD
{
    public partial class VoronoiForm : Form
    {
        private static readonly int POINT_CLOUD_SIZE = 420;
        private static readonly int NUMBER_OF_ITERATIONS = 100;

        private readonly Graphics graphicsContext;
        private readonly PictureBox canvasBox = new();


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
            //Point3D basePoint = new(3,3,3);
            //Point3D normalVector = new(1,1,1);
            //DrawPlane(normalVector, basePoint);



            List<Point3D> randomPointCloud = Point3DGenerator.GenerateRandomPoints(POINT_CLOUD_SIZE, WIDTH, HEIGHT);
            Point3D basePoint = new(0,0,0);
            Point3D normalVector = new(1,1,1);
            Plane2D plane2D = new(normalVector, basePoint);
            PointTransformation pointTransformation = CreatePointTransformation(plane2D);

            List<VoronoiPoint2D> projectedPointloud = ProjectPointCloudOnGround(plane2D, randomPointCloud, pointTransformation);
            DelaunayTriangulation triangulationProcess = new DelaunayTriangulation();
            ISet<DelaunayTriangle> triangulation = triangulationProcess.CreateTriangulation(projectedPointloud);
            VoronoiDiagram voronoi = new();
            Dictionary<VoronoiPoint2D, VoronoiCell3D> diagram = voronoi.CreateVoronoiDiagram(triangulation);
            diagram = CenterVoronoiDiagram(diagram, NUMBER_OF_ITERATIONS);
            List<VoronoiCell3D> projectedDiagram = ProjectVoronoiCellsReversively(diagram, pointTransformation);


            Point3D screenCenter = new(WIDTH / 2, HEIGHT / 2, 0);
            
            //Debug.WriteLine(averagePoint.ToString());
            
            projectedDiagram = ProjectVoronoiCellsOnPlane(basePoint, normalVector, projectedDiagram);
            Point3D averagePoint = CalculateAveragePoint(projectedDiagram);
            Point3D translationVector = new(screenCenter.X - averagePoint.X, screenCenter.Y - averagePoint.Y, 0);
            Translation transform = new();
            transform.SetTranslation(translationVector);
            projectedDiagram = ProjectVoronoiCells(projectedDiagram, transform);
            DrawVoronoiCells(projectedDiagram);
        }

        private List<VoronoiCell3D> ProjectVoronoiCellsOnPlane(Point3D basePoint, Point3D normalVector, List<VoronoiCell3D> voronoiCells)
        {
            float angleX = 60;
            float angleY = 0;
            float angleZ = 0;
            List<VoronoiCell3D> projectedCells = new();
            Func<Point3D, Point3D> rotate = (p) =>
            {
                // Convert angles to radians
                float radX = angleX * (float)Math.PI / 180;
                float radY = angleY * (float)Math.PI / 180;
                float radZ = angleZ * (float)Math.PI / 180;

                // Rotate around X axis
                double y1 = p.Y * (float)Math.Cos(radX) - p.Z * (float)Math.Sin(radX);
                double z1 = p.Y * (float)Math.Sin(radX) + p.Z * (float)Math.Cos(radX);

                // Rotate around Y axis
                double x2 = p.X * (float)Math.Cos(radY) + z1 * (float)Math.Sin(radY);
                double z2 = -p.X * (float)Math.Sin(radY) + z1 * (float)Math.Cos(radY);

                // Rotate around Z axis
                double x3 = x2 * (float)Math.Cos(radZ) - y1 * (float)Math.Sin(radZ);
                double y3 = x2 * (float)Math.Sin(radZ) + y1 * (float)Math.Cos(radZ);

                return new Point3D(x3, y3, z2);
            };

            Func<Point3D, Point> project = (p) =>
            {
                float scale = 1; // Scale factor for better visualization
                return new Point((int)(WIDTH / 2 + p.X * scale),
                    (int)(HEIGHT / 2 - p.Y * scale)
                );
            };
            foreach (VoronoiCell3D voronoiCell in voronoiCells)
            {
                VoronoiCell3D projectedCell = ProjectVoronoiCell(basePoint, normalVector, rotate, project, voronoiCell); 
                projectedCells.Add(projectedCell);
            }

            DrawAxis(project, rotate, Color.Red, new Point3D(-500, 0, 0), new Point3D(500, 0, 0));
            DrawAxis(project, rotate, Color.Green, new Point3D(0, -500, 0), new Point3D(0, 500, 0));
            DrawAxis(project, rotate, Color.Blue, new Point3D(0, -500, -500), new Point3D(0, -500, -500));


            return projectedCells;
        }

        private VoronoiCell3D ProjectVoronoiCell(Point3D basePoint, Point3D normalVector, Func<Point3D, Point3D> rotate, Func<Point3D, Point> project, VoronoiCell3D voronoiCell)
        {
            Point3D center = voronoiCell.getCenter();
            List<Edge3D> edges = voronoiCell.GetEdges();

            List<Edge3D> projectedEdges = new();
            foreach (Edge3D edge in edges)
            {
                Edge3D projectedEdge = ProjectVoronoiEdge(basePoint, normalVector, rotate, project, edge);
                projectedEdges.Add(projectedEdge);
            }
            Point projectCenter = ProjectPoint(basePoint, normalVector, rotate, project, center);
            Point3D cntr = new(projectCenter.X, projectCenter.Y, 0);

            return new(cntr, projectedEdges);
        }

        private Point ProjectPoint(Point3D basePoint, Point3D normalVector, Func<Point3D, Point3D> rotate, Func<Point3D, Point> project, Point3D pointToProject)
        {
            double z = (normalVector.X * (pointToProject.X - basePoint.X) + normalVector.Y * (pointToProject.Y - basePoint.Y)) / -normalVector.Z + basePoint.Z;
            Point3D p = new Point3D(pointToProject.X, pointToProject.Y, z);
            Point3D rotatedPoint = rotate(p);
            Debug.WriteLine(project(rotatedPoint));
            return project(rotatedPoint);
        }

        private Edge3D ProjectVoronoiEdge(Point3D basePoint, Point3D normalVector, Func<Point3D, Point3D> rotate, Func<Point3D, Point> project, Edge3D edgeToProject)
        {
            Point3D start = edgeToProject.startingPoint;
            Point3D end = edgeToProject.endingPoint;
            
            Point projectedStart = ProjectPoint(basePoint, normalVector, rotate, project, start);  
            Point projectedEnd = ProjectPoint(basePoint, normalVector, rotate, project, end);

            Point3D pS = new(projectedStart.X, projectedStart.Y, 0);
            Point3D pE = new(projectedEnd.X, projectedEnd.Y, 0);
            return new(pS, pE);
        }

        private void DrawPlane(Point3D normal, Point3D point)
        {
            float angleX = 0;
            float angleY = 45;
            float angleZ = 0;

            Plane2D plane = new(normal, point);
            Function planeEquation = plane.equation;

            Func<Point3D, Point3D> rotate = (p) =>
            {
                // Convert angles to radians
                float radX = angleX * (float)Math.PI / 180;
                float radY = angleY * (float)Math.PI / 180;
                float radZ = angleZ * (float)Math.PI / 180;

                // Rotate around X axis
                double y1 = p.Y * (float)Math.Cos(radX) - p.Z * (float)Math.Sin(radX);
                double z1 = p.Y * (float)Math.Sin(radX) + p.Z * (float)Math.Cos(radX);

                // Rotate around Y axis
                double x2 = p.X * (float)Math.Cos(radY) + z1 * (float)Math.Sin(radY);
                double z2 = -p.X * (float)Math.Sin(radY) + z1 * (float)Math.Cos(radY);

                // Rotate around Z axis
                double x3 = x2 * (float)Math.Cos(radZ) - y1 * (float)Math.Sin(radZ);
                double y3 = x2 * (float)Math.Sin(radZ) + y1 * (float)Math.Cos(radZ);

                return new Point3D(x3, y3, z2);
            };

            Func<Point3D, Point> project = (p) =>
            {
                float scale = 30; // Scale factor for better visualization
                return new Point(
                    (int)(WIDTH / 2 + p.X * scale),
                    (int)(HEIGHT / 2 - p.Y * scale)
                );
            };

            // Draw axes
            DrawAxis(project, rotate, Color.Red, new Point3D(-5, 0, 0), new Point3D(5, 0, 0));
            DrawAxis(project, rotate, Color.Green, new Point3D(0, -5, 0), new Point3D(0, 5, 0));
            DrawAxis(project, rotate, Color.Blue, new Point3D(0, -5, -5), new Point3D(0, -5, -5));

            // Calculate points on the plane
            for (float x = -5; x <= 5; x += 0.1f)
            {
                for (float y = -5; y <= 5; y += 0.1f)
                {
                    double z = (normal.X * (x - point.X) + normal.Y * (y - point.Y)) / -normal.Z + point.Z;
                    Point3D p = new Point3D(x, y, z);
                    Point3D rotatedPoint = rotate(p);
                    Point p2D = project(rotatedPoint);
                    if (p2D.X >= 0 && p2D.X < WIDTH && p2D.Y >= 0 && p2D.Y < HEIGHT)
                    {
                        ((Bitmap)canvasBox.Image).SetPixel(p2D.X, p2D.Y, Color.Black);
                    }
                }
            }
        }

        private void DrawAxis(Func<Point3D, Point> project, Func<Point3D, Point3D> rotate, Color color, Point3D start, Point3D end)
        {
            using (Pen pen = new Pen(color, 2))
            {
                Point3D rotatedStart = rotate(start);
                Point3D rotatedEnd = rotate(end);
                graphicsContext.DrawLine(pen, project(rotatedStart), project(rotatedEnd));

                for (float i = -5; i <= 5; i++)
                {
                    Point3D tickStart = new Point3D(start.X + i * (end.X - start.X) / 10, start.Y + i * (end.Y - start.Y) / 10, start.Z + i * (end.Z - start.Z) / 10);
                    Point3D tickEnd = new Point3D(tickStart.X, tickStart.Y, tickStart.Z);
                    tickEnd.X += 0.1f * (end.Y - start.Y);
                    tickEnd.Y += 0.1f * (end.Z - start.Z);
                    tickEnd.Z += 0.1f * (end.X - start.X);
                    Point3D rotatedTickStart = rotate(tickStart);
                    Point3D rotatedTickEnd = rotate(tickEnd);
                    //graphicsContext.DrawLine(pen, project(rotatedTickStart), project(rotatedTickEnd));
                }
            }
        }

        private Point3D CalculateAveragePoint(List<VoronoiCell3D> cells)
        {
            double totalX = 0;
            double totalY = 0;
            double totalZ = 0;
            int bound = Math.Max(WIDTH, HEIGHT) * 2;
            foreach (VoronoiCell3D cell in cells) 
            {
                Point3D center = cell.getCenter();
                if (Math.Abs(center.X) > bound || Math.Abs(center.Y) > bound || Math.Abs(center.Z) > bound)
                {
                    continue;
                }
                totalX += center.X;
                totalY += center.Y;
                totalZ += center.Z;
            }

            double finalX = totalX / cells.Count;
            double finalY = totalY / cells.Count;
            double finalZ = totalZ / cells.Count;

            return new(finalX, finalY, finalZ);
        }

        private PointTransformation CreatePointTransformation(Plane2D plane)
        {
            Point3D groundNormalVector = new(0,0,1);
            Point3D initialPoint = new(0, 0, 0);
            Point3D normalVector = plane.normalVector;
            Plane2D groundPlane = new(groundNormalVector, initialPoint);

            Point3D rotationAxisVector = normalVector.CalculateCrossProduct(groundNormalVector);
            double angle = plane.CalculateAngle(groundPlane);
            
            Rotation rotation = new(rotationAxisVector, (float) angle);
            Translation translation = new Translation();
            Point3D? intersectPoint = CalculateIntersectVectorWithGround(plane);


            if (intersectPoint != null)
            {
                translation.SetTranslation(intersectPoint.Multiply(-1));
            }

            return new(translation, rotation);

        }

        private List<VoronoiPoint2D> ProjectPointCloudOnGround(Plane2D plane, List<Point3D> pointcloud, PointTransformation pointTransformation)
        {
            List<VoronoiPoint2D> projectedPointcloud = new();
            double epsilon = 1e-4;
            foreach (Point3D point in pointcloud)
            {
                Point3D projectedPointOnPlane = ProjectPointOnPlane(plane, point);
                Point3D projectedPoint = pointTransformation.Transform(projectedPointOnPlane);
               
                if (projectedPoint.Z < epsilon)
                {
                    VoronoiPoint2D projectedPointWithoutZAxe = new(projectedPoint.X, projectedPoint.Y);
                    projectedPointcloud.Add(projectedPointWithoutZAxe);
                } else
                {
                    MessageBox.Show(projectedPoint.ToString());
                    Debug.WriteLine("Point projection did not go well - point " + point.ToString() + " was not displayed on ground properly.");
                }
            }

            return projectedPointcloud;
        }

        private Point3D ProjectPointOnPlane(Plane2D plane, Point3D point)
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

            return projectedPoint;
        }
        private Point3D? CalculateIntersectVectorWithGround(Plane2D plane)
        {
            double a = plane.equation.A;
            double b = plane.equation.B;
            double d = plane.equation.D;
            if (a == 0 && b == 0)
            {
                return null;
            }
            if (b == 0)
            {
                return new(-d / a, 0, 0);
            }
            double x = 2;
            double y = (-d - a*x) / b;
            return new(x, y, 0);
        }

        private void Draw3DPoint(Brush brush, Point3D point)
        {
            int pointWidth = 4;
            int pointHeight = 4;
            if (point.X > 0 && point.Y > 0)
            {
                graphicsContext.FillEllipse(brush, (int)point.X - pointWidth / 2, (int)point.Y - pointHeight / 2, pointWidth, pointHeight);
            }
        }

        private Dictionary<VoronoiPoint2D, VoronoiCell3D> CalculateVoronoiDiagram(List<VoronoiPoint2D> points)
        {
            DelaunayTriangulation delaunayTriangulation = new();
            ISet<DelaunayTriangle> triangulation = delaunayTriangulation.CreateTriangulation(points);
            
            VoronoiDiagram voronoi = new();
            return voronoi.CreateVoronoiDiagram(triangulation);
        }

        private Dictionary<VoronoiPoint2D, VoronoiCell3D> CenterVoronoiDiagram(Dictionary<VoronoiPoint2D, VoronoiCell3D> voronoiDiagram, int numberOfIterations)
        {
            for (int i = 0; i < numberOfIterations; i++)
            {
                List<VoronoiPoint2D> newPoints = new();
                foreach (VoronoiCell3D cell in voronoiDiagram.Values)
                {
                    newPoints.Add(CalculateCenterOfGravity(cell));
                }
                voronoiDiagram = CalculateVoronoiDiagram(newPoints);
            }

            return voronoiDiagram;
        }

        private List<VoronoiCell3D> ProjectVoronoiCellsReversively(Dictionary<VoronoiPoint2D, VoronoiCell3D> diagram, RevertableTransform transform)
        {
            List<VoronoiCell3D> projectedCells = new();
            foreach (VoronoiCell3D voronoiCell in diagram.Values)
            {
                Point3D centerPoint = voronoiCell.getCenter();
                List<Edge3D> cellEdges = voronoiCell.GetEdges();

                Point3D projectedCenter = transform.RevertTransform(centerPoint);
                VoronoiCell3D projectedCell = new(projectedCenter);

                foreach (Edge3D cellEdge in cellEdges)
                {
                    Point3D projectedStartingVertex = transform.RevertTransform(cellEdge.startingPoint);
                    Point3D projectedEndingVertex = transform.RevertTransform(cellEdge.endingPoint);
                    Edge3D projectedEdge = new(projectedStartingVertex, projectedEndingVertex);
                    projectedCell.AddEdge(projectedEdge);
                }

                projectedCells.Add(projectedCell);
            }

            return projectedCells;
        }

        private List<VoronoiCell3D> ProjectVoronoiCells(List<VoronoiCell3D> cells, Transform transform)
        {
            List<VoronoiCell3D> projectedCells = new();
            foreach (VoronoiCell3D voronoiCell in cells)
            {
                Point3D centerPoint = voronoiCell.getCenter();
                List<Edge3D> cellEdges = voronoiCell.GetEdges();

                Point3D projectedCenter = transform.Transform(centerPoint);
                VoronoiCell3D projectedCell = new(projectedCenter);

                foreach (Edge3D cellEdge in cellEdges)
                {
                    Point3D projectedStartingVertex = transform.Transform(cellEdge.startingPoint);
                    Point3D projectedEndingVertex = transform.Transform(cellEdge.endingPoint);
                    Edge3D projectedEdge = new(projectedStartingVertex, projectedEndingVertex);
                    projectedCell.AddEdge(projectedEdge);
                }

                projectedCells.Add(projectedCell);
            }

            return projectedCells;
        }

        private VoronoiPoint2D CalculateCenterOfGravity(VoronoiCell3D voronoiCell)
        {
            ISet<Point3D> vertices = voronoiCell.getVertices();
            int verticesSize = vertices.Count;
            double findingX = 0;
            double findingY = 0;
            foreach (Point3D vertex in vertices)
            {
                if (Math.Abs(vertex.X) > WIDTH || Math.Abs(vertex.Y) > HEIGHT)
                {
                    Point3D centerPoint = voronoiCell.getCenter();
                    return new(centerPoint.X, centerPoint.Y);
                }

                findingX += vertex.X;
                findingY += vertex.Y;
            }

            return new(findingX / verticesSize, findingY / verticesSize);
        }

        private void DrawVoronoiCells(ICollection<VoronoiCell3D> voronoiCells)
        {
            Pen pen = new(Color.Blue);
            Brush pointPen = new SolidBrush(Color.Black);
            foreach (VoronoiCell3D cell in voronoiCells)
            {             
                DrawVoronoiCell(pen, pointPen, cell);
            }
        }

        private void DrawVoronoiCell(Pen pen, Brush pointPen, VoronoiCell3D cell)
        {
            foreach (Edge3D vertex in cell.GetEdges())
            {
                graphicsContext.DrawLine(pen,
                    new((int) vertex.startingPoint.X, (int) vertex.startingPoint.Y),
                    new((int) vertex.endingPoint.X, (int) vertex.endingPoint.Y));
            }

            Draw3DPoint(pointPen, cell.getCenter());
        }

        private void DrawTriangles(ICollection<DelaunayTriangle> triangles)
        {
            Pen pen = new Pen(Color.Red);
            foreach(DelaunayTriangle triangle in triangles)
            {
                DrawTriangle(pen, triangle);
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
