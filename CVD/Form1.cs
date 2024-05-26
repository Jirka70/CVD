using System.Diagnostics;
using System.Numerics;
using CVD.delaunay;
using CVD.shape;
using CVD.transform;
using CVD.transforms;
using CVD.util;
using CVD.voronoi;

namespace CVD
{
    public partial class VoronoiForm : Form
    {
        private static readonly int POINT_CLOUD_SIZE = 420;
        private static readonly int NUMBER_OF_ITERATIONS = 10;

        private readonly TextBox xAngleTextBox = new();
        private readonly TextBox yAngleTextBox = new();
        private readonly TextBox zAngleTextBox = new();
        private readonly TextBox scaleBox = new();
        private readonly Button applyButton = new();

        private readonly Graphics graphicsContext;
        private readonly PictureBox canvasBox = new();
        private readonly List<VoronoiCell3D> diagram;


        public VoronoiForm()
        {
            InitializeComponent();
            canvasBox.Image = new Bitmap(WIDTH, HEIGHT);
            canvasBox.Size = new Size(WIDTH, HEIGHT);
            canvasBox.Location = new Point(0, 0);
            Controls.Add(canvasBox);

            PictureBox box = new PictureBox();
            box.Size = new Size(100, 100);
            box.Location = new Point(0, 0);
            box.Image = new Bitmap(100, 100);
            Controls.Add(box);
            Bitmap image = (Bitmap)canvasBox.Image;
            graphicsContext = Graphics.FromImage(image);

            Label xAngleLabel = new();
            Label yAngleLabel = new();
            Label zAngleLabel = new();
            Label scaleLabel = new();

            SetUpLabel(xAngleLabel, "X Angle", 0);
            xAngleTextBox.Location = new Point(200, 200);
            Controls.Add(xAngleTextBox);

            SetUpLabel(xAngleLabel, "X Angle", 0);
            SetUpTextBox(xAngleTextBox, 0);

            SetUpLabel(yAngleLabel, "Y Angle", 1);
            SetUpTextBox(yAngleTextBox, 1);

            SetUpLabel(zAngleLabel, "Z Angle", 2);
            SetUpTextBox(zAngleTextBox, 2);

            SetUpLabel(scaleLabel, "Scale", 3);
            SetUpTextBox(scaleBox, 3);

            SetUpButton(applyButton, "Apply", 4);
            Point3D basePoint = new(100,100,100);
            Point3D normalVector = new(1, 1, 1);
            Plane2D mainPlane = new(normalVector, basePoint);
            

            // Set initial angles and hook up the event
            xAngleTextBox.Text = "0";
            yAngleTextBox.Text = "0";
            zAngleTextBox.Text = "0";
            scaleBox.Text = "0.5";

            List<Point3D> randomPointCloud = Point3DGenerator.GenerateRandomPoints(POINT_CLOUD_SIZE, WIDTH, HEIGHT); 


            diagram = CalculateVoronoiDiagram(randomPointCloud, mainPlane, normalVector, basePoint);
            List<VoronoiCell3D> projectedDiagram = ApplyRotation(diagram, 0,0,0,.5f, mainPlane, normalVector, basePoint);
            DrawVoronoiCells(projectedDiagram);
            applyButton.Click += (o, e) => ApplyButton_Click(o, e, mainPlane, normalVector, basePoint, diagram);
            //VoronoiProcess();
        }

        private void SetUpLabel(Label label, string labelText, int rowNumber)
        {
            int rowHeight = 30;
            label.Text = labelText;
            label.Location = new Point(10, rowHeight * rowNumber + 10);
            label.AutoSize = true;
            label.BackColor = Color.Transparent; // Make the label background transparent
            canvasBox.Controls.Add(label);
        }

        private void ApplyButton_Click(object? sender, EventArgs e, Plane2D mainPlane, Point3D normalVector, Point3D basePoint, List<VoronoiCell3D> voronoiCells)
        {
            if (float.TryParse(xAngleTextBox.Text, out float angleX) &&
                float.TryParse(yAngleTextBox.Text, out float angleY) &&
                float.TryParse(zAngleTextBox.Text, out float angleZ) &&
                float.TryParse(scaleBox.Text, out float scale))
            {

                graphicsContext.Clear(canvasBox.BackColor);
               
                List<VoronoiCell3D> rotatedCells = ApplyRotation(voronoiCells, angleX, angleY, angleZ, scale, mainPlane, normalVector, basePoint);
                DrawVoronoiCells(rotatedCells);
                canvasBox.Refresh();
            }
            else
            {
                MessageBox.Show("Please enter valid angle values.", "Invalid Input", MessageBoxButtons.OK, MessageBoxIcon.Error);
            }
        }

        private List<VoronoiCell3D> ApplyRotation(List<VoronoiCell3D> voronoiCells, float angleX, float angleY, float angleZ, float scale, Plane2D mainPlane, Point3D normalVector, Point3D basePoint)
        {
            Func<Point3D, Point3D> rotate = (p) =>
            {
                float radX = angleX * (float)Math.PI / 180;
                float radY = angleY * (float)Math.PI / 180;
                float radZ = angleZ * (float)Math.PI / 180;

                double y1 = p.Y * (float)Math.Cos(radX) - p.Z * (float)Math.Sin(radX);
                double z1 = p.Y * (float)Math.Sin(radX) + p.Z * (float)Math.Cos(radX);

                double x2 = p.X * (float)Math.Cos(radY) + z1 * (float)Math.Sin(radY);
                double z2 = -p.X * (float)Math.Sin(radY) + z1 * (float)Math.Cos(radY);

                double x3 = x2 * (float)Math.Cos(radZ) - y1 * (float)Math.Sin(radZ);
                double y3 = x2 * (float)Math.Sin(radZ) + y1 * (float)Math.Cos(radZ);

                return new Point3D(x3, y3, z2);
            };

            Func<Point3D, Point3D> project = (p) =>
            {
                double x = WIDTH / 2 + p.X * scale;
                double y = HEIGHT / 2 - p.Y * scale;
                return new(x, y, 0);
            };

            return ProjectVoronoiCellsOnPlane(mainPlane, basePoint, rotate, project, normalVector, voronoiCells);
        }


        private void SetUpTextBox(TextBox textBox, int rowNumber)
        {
            int rowHeight = 30;
            textBox.Location = new Point(80, rowHeight * rowNumber + 10);
            textBox.Size = new Size(100, 20);
            canvasBox.Controls.Add(textBox);
        }

        private void SetUpButton(Button button, string buttonText, int rowNumber)
        {
            int rowHeight = 30;
            button.Text = buttonText;
            button.Location = new Point(10, rowHeight * rowNumber + 10);
            button.Size = new Size(170, 30);
            canvasBox.Controls.Add(button);
        }


        private List<VoronoiCell3D> CalculateVoronoiDiagram(List<Point3D> pointCloud, Plane2D mainPlane, Point3D normalVector, Point3D basePoint)
        {
            PointTransformation pointTransformation = CreatePointTransformation(mainPlane);

            List<VoronoiPoint2D> projectedPointloud = ProjectPointCloudOnGround(mainPlane, pointCloud, pointTransformation);
            DelaunayTriangulation triangulationProcess = new DelaunayTriangulation();
            ISet<DelaunayTriangle> triangulation = triangulationProcess.CreateTriangulation(projectedPointloud);
            VoronoiDiagram voronoi = new();
            Dictionary<VoronoiPoint2D, VoronoiCell3D> diagram = voronoi.CreateVoronoiDiagram(triangulation);
            diagram = CenterVoronoiDiagram(diagram, NUMBER_OF_ITERATIONS);
            List<VoronoiCell3D> projectedDiagram = ProjectVoronoiCellsReversively(diagram, pointTransformation);

            return projectedDiagram;
        }

        private List<VoronoiCell3D> ProjectVoronoiCellsOnPlane(Plane2D plane, Point3D basePoint, Func<Point3D, Point3D> rotate, Func<Point3D, Point3D> project, Point3D normalVector, List<VoronoiCell3D> voronoiCells)
        {
            Point3D[] corners = CalculatePlaneCorners(basePoint, normalVector, 10f, 10f);
            Function planeEquation = plane.equation;
            
            List<VoronoiCell3D> projectedCells = new();
            
            foreach (VoronoiCell3D voronoiCell in voronoiCells)
            {
                VoronoiCell3D projectedCell = ProjectVoronoiCell(basePoint, normalVector, corners, rotate, project, voronoiCell); 
                projectedCells.Add(projectedCell);
            }

            DrawAxis(project, rotate, Color.Red, new Point3D(-500, 0, 0), new Point3D(500, 0, 0), "X");
            DrawAxis(project, rotate, Color.Green, new Point3D(0, -500, 0), new Point3D(0, 500, 0), "Y");
            DrawAxis(project, rotate, Color.Blue, new Point3D(0, 0, -500), new Point3D(0, 0, 500), "Z");
            //DrawPlane(normalVector, basePoint, rotate, project);

            return projectedCells;
        }

        private VoronoiCell3D ProjectVoronoiCell(Point3D basePoint, Point3D normalVector, Point3D[] corners, Func<Point3D, Point3D> rotate, Func<Point3D, Point3D> project, VoronoiCell3D voronoiCell)
        {
            Point3D center = voronoiCell.getCenter();
            List<Edge3D> edges = voronoiCell.GetEdges();

            List<Edge3D> projectedEdges = new();
            foreach (Edge3D edge in edges)
            {
                Edge3D projectedEdge = ProjectVoronoiEdge(basePoint, normalVector, corners, rotate, project, edge);
                projectedEdges.Add(projectedEdge);
            }

            Point3D projectedCenter = ProjectPoint(basePoint, normalVector, corners, rotate, project, center);

            return new(projectedCenter, projectedEdges);
        }

        private Point3D ProjectPoint(Point3D basePoint, Point3D normalVector, Point3D[] corners, Func<Point3D, Point3D> rotate, Func<Point3D, Point3D> project, Point3D pointToProject)
        {
            double z = (normalVector.X * (pointToProject.X - basePoint.X) + normalVector.Y * (pointToProject.Y - basePoint.Y)) / -normalVector.Z + basePoint.Z;
            Point3D p = new(pointToProject.X, pointToProject.Y, z);
            Point3D projectedPoint = project(rotate(p));
            return projectedPoint;
        }

        private Edge3D ProjectVoronoiEdge(Point3D basePoint, Point3D normalVector, Point3D[] corners, Func<Point3D, Point3D> rotate, Func<Point3D, Point3D> project, Edge3D edgeToProject)
        {
            Point3D start = edgeToProject.startingPoint;
            Point3D end = edgeToProject.endingPoint;


            Point3D projectedStart = ProjectPoint(basePoint, normalVector, corners, rotate, project, start);
            Point3D projectedEnd = ProjectPoint(basePoint, normalVector, corners, rotate, project, end);

            return new(projectedStart, projectedEnd);
        }

        private void DrawPlane(Point3D normal, Point3D point, Func<Point3D, Point3D> rotate, Func<Point3D, Point3D> project)
        {
            Point3D[] corners = CalculatePlaneCorners(point, normal, WIDTH, HEIGHT);

            Plane2D plane = new(normal, point);

            

            Color semiTransparentColor = Color.FromArgb(50, Color.LightBlue);
            Brush semiTransparentBrush = new SolidBrush(semiTransparentColor);

            Point3D v1 = new(corners[1].X - corners[0].X, corners[1].Y - corners[0].Y, corners[1].Z - corners[0].Z);
            Point3D v2 = new(corners[3].X - corners[0].X, corners[3].Y - corners[0].Y, corners[3].Z - corners[0].Z);

            for (float i = 0; i <= 1; i += 0.005f)
            {
                for (float j = 0; j <= 1; j += 0.005f)
                {
                    Point3D p = new(
                        corners[0].X + i * v1.X + j * v2.X,
                        corners[0].Y + i * v1.Y + j * v2.Y,
                        corners[0].Z + i * v1.Z + j * v2.Z
                    );

                    Point3D rotatedPoint = rotate(p);
                    Point3D p2D = project(rotatedPoint);

                    if (p2D.X >= 0 && p2D.X < WIDTH && p2D.Y >= 0 && p2D.Y < HEIGHT)
                    {
                        graphicsContext.FillRectangle(semiTransparentBrush, (int) p2D.X, (int) p2D.Y, 2, 2); // Draw point with semi-transparent color
                    }
                }
            }
        }
        private Point3D[] CalculatePlaneCorners(Point3D center, Point3D normal, float width, float height)
        {
            // Calculate two vectors in the plane
            Point3D v1, v2;
            if (normal.X != 0 || normal.Y != 0)
            {
                v1 = new(-normal.Y, normal.X, 0);
            }
            else
            {
                v1 = new(0,-normal.Z, normal.Y);
            }
            v2 = CrossProduct(normal, v1);

            // Normalize and scale vectors
            v1 = Normalize(v1);
            v2 = Normalize(v2);
            v1 = Scale(v1, width / 2);
            v2 = Scale(v2, height / 2);

            // Calculate corners
            return new Point3D[]
            {
            new(center.X + v1.X + v2.X, center.Y + v1.Y + v2.Y, center.Z + v1.Z + v2.Z),
            new(center.X + v1.X - v2.X, center.Y + v1.Y - v2.Y, center.Z + v1.Z - v2.Z),
            new(center.X - v1.X - v2.X, center.Y - v1.Y - v2.Y, center.Z - v1.Z - v2.Z),
            new(center.X - v1.X + v2.X, center.Y - v1.Y + v2.Y, center.Z - v1.Z + v2.Z)
            };
        }

        private Point3D CrossProduct(Point3D a, Point3D b)
        {
            return new(
                a.Y * b.Z - a.Z * b.Y,
                a.Z * b.X - a.X * b.Z,
                a.X * b.Y - a.Y * b.X
            );
        }

        private Point3D Normalize(Point3D v)
        {
            float length = (float)Math.Sqrt(v.X * v.X + v.Y * v.Y + v.Z * v.Z);
            return new(v.X / length, v.Y / length, v.Z / length);
        }

        private Point3D Scale(Point3D v, float scale)
        {
            return new(v.X * scale, v.Y * scale, v.Z * scale);
        }

        private void DrawAxis(Func<Point3D, Point3D> project, Func<Point3D, Point3D> rotate, Color color, Point3D start, Point3D end, string axisName)
        {
            Point3D startProjected = project(rotate(start));
            Point3D endProjected = project(rotate(end));

            Point s = new((int)startProjected.X, (int)startProjected.Y);
            Point e = new((int)endProjected.X, (int)endProjected.Y);
            graphicsContext.DrawLine(new Pen(color, 2), s,e);

            DrawAxisLabel(project, rotate, color, start, axisName);
            DrawAxisLabel(project, rotate, color, end, axisName);

            int interval = 100; // Define interval for the labels
            for (int i = (int)start.X+interval; i < end.X - interval; i += interval)
            {
                if (axisName == "X")
                    DrawAxisLabel(project, rotate, color, new Point3D(i, start.Y, start.Z), i.ToString());
            }
            for (int i = (int)start.Y + interval; i < end.Y - interval; i += interval)
            {
                if (axisName == "Y")
                    DrawAxisLabel(project, rotate, color, new Point3D(start.X, i, start.Z), i.ToString());
            }
            for (int i = (int)start.Z + interval; i <= end.Z - interval; i += interval)
            {
                if (axisName == "Z")
                    DrawAxisLabel(project, rotate, color, new Point3D(start.X, start.Y, i), i.ToString());
            }
        }

        private void DrawAxisLabel(Func<Point3D, Point3D> project, Func<Point3D, Point3D> rotate, Color color, Point3D position, string labelValue)
        {
            Point3D projectedPosition = project(rotate(position));
            Debug.WriteLine(projectedPosition.X + " " + projectedPosition.Y);
                
            Point pos = new((int)projectedPosition.X, (int)projectedPosition.Y);

            graphicsContext.DrawString(labelValue, new Font("Arial", 8), new SolidBrush(color), pos);
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


        private VoronoiPoint2D CalculateCenterOfGravity(VoronoiCell3D voronoiCell)
        {
            ISet<Point3D> vertices = voronoiCell.getVertices();
            int verticesSize = vertices.Count;
            double findingX = 0;
            double findingY = 0;
            foreach (Point3D vertex in vertices)
            {
                if (Math.Abs(vertex.X) > 10000 || Math.Abs(vertex.Y) > 10000)
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
            Pen pen = new(Color.Brown);
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
    }
}
