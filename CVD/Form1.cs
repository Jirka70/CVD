using System.Diagnostics;
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
        private static readonly int NUMBER_OF_ITERATIONS = 100;

        private readonly TextBox xAngleTextBox = new();
        private readonly TextBox yAngleTextBox = new();
        private readonly TextBox zAngleTextBox = new();
        private readonly Button applyButton = new();

        private readonly Graphics graphicsContext;
        private readonly PictureBox canvasBox = new();
        private readonly Dictionary<VoronoiPoint2D, VoronoiCell3D> diagram;


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

            SetUpLabel(xAngleLabel, "X Angle", 0);
            xAngleTextBox.Location = new Point(200, 200);
            xAngleTextBox.Text = "piko";
            Controls.Add(xAngleTextBox);

            SetUpLabel(xAngleLabel, "X Angle", 0);
            SetUpTextBox(xAngleTextBox, 0);

            SetUpLabel(yAngleLabel, "Y Angle", 1);
            SetUpTextBox(yAngleTextBox, 1);

            SetUpLabel(zAngleLabel, "Z Angle", 2);
            SetUpTextBox(zAngleTextBox, 2);

            SetUpButton(applyButton, "Apply", 3);

            applyButton.Click += ApplyButton_Click;



            Text = "3D Voronoi Projection";

            // Set initial angles and hook up the event
            xAngleTextBox.Text = "30";
            yAngleTextBox.Text = "30";
            zAngleTextBox.Text = "30";

            VoronoiProcess();
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

        private void ApplyButton_Click(object sender, EventArgs e)
        {
            if (float.TryParse(xAngleTextBox.Text, out float x) &&
                float.TryParse(yAngleTextBox.Text, out float y) &&
                float.TryParse(zAngleTextBox.Text, out float z))
            {
                //RedrawGraphics(x, y, z);
            }
            else
            {
                MessageBox.Show("Please enter valid angle values.", "Invalid Input", MessageBoxButtons.OK, MessageBoxIcon.Error);
            }
        }

        private void ApplyRotation(float angleX, float angleY, float angleZ)
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
                float scale = .5f;
                double x = WIDTH / 2 + p.X * scale;
                double y = HEIGHT / 2 - p.Y * scale;
                return new(
                    x,
                    y, 0);
            };

            graphicsContext.Clear(Color.White);

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

        private void VoronoiProcess()
        {

            Point3D basePoint = new(1, 1, 1);
            Point3D normalVector = new(1, 1, 1);
            Plane2D plane2D = new(normalVector, basePoint);

            Function planeEquation = plane2D.equation;
            float angleX = 0;
            float angleY = 30;
            float angleZ = 0;
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
                float scale = .5f;
                double x = WIDTH / 2 + p.X * scale;
                double y = HEIGHT / 2 - p.Y * scale;
                return new(
                    x,
                    y, 0);
            };

            List<Point3D> randomPointCloud = Point3DGenerator.GenerateRandomPoints(POINT_CLOUD_SIZE, WIDTH, HEIGHT);
           
            Debug.WriteLine(plane2D.equation);
            PointTransformation pointTransformation = CreatePointTransformation(plane2D);

            List<VoronoiPoint2D> projectedPointloud = ProjectPointCloudOnGround(plane2D, randomPointCloud, pointTransformation);
            DelaunayTriangulation triangulationProcess = new DelaunayTriangulation();
            ISet<DelaunayTriangle> triangulation = triangulationProcess.CreateTriangulation(projectedPointloud);
            VoronoiDiagram voronoi = new();
            Dictionary<VoronoiPoint2D, VoronoiCell3D> diagram = voronoi.CreateVoronoiDiagram(triangulation);
            diagram = CenterVoronoiDiagram(diagram, NUMBER_OF_ITERATIONS);
            List<VoronoiCell3D> projectedDiagram = ProjectVoronoiCellsReversively(diagram, pointTransformation);


            projectedDiagram = ProjectVoronoiCellsOnPlane(plane2D, basePoint, rotate, project, normalVector, projectedDiagram);
            Point3D averagePoint = CalculateAveragePoint(projectedDiagram);
            
            DrawVoronoiCells(projectedDiagram);
            DrawPlane(normalVector, basePoint, rotate, project);
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

            DrawAxis(project, rotate, Color.Red, new Point3D(-500, 0, 0), new Point3D(500, 0, 0));
            DrawAxis(project, rotate, Color.Green, new Point3D(0, -500, 0), new Point3D(0, 500, 0));
            DrawAxis(project, rotate, Color.Blue, new Point3D(0, -500, -500), new Point3D(0, -500, -500));


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
//            Point3D cntr = new(projectedCenter.X, projectedCenter.Y, 0);


            return new(projectedCenter, projectedEdges);



        }

        private Point3D ProjectPoint(Point3D basePoint, Point3D normalVector, Point3D[] corners, Func<Point3D, Point3D> rotate, Func<Point3D, Point3D> project, Point3D pointToProject)
        {
            double z = (normalVector.X * (pointToProject.X - basePoint.X) + normalVector.Y * (pointToProject.Y - basePoint.Y)) / -normalVector.Z + basePoint.Z;
            Point3D p = new(pointToProject.X, pointToProject.Y, z);
            Point3D projectedPoint = project(rotate(p));
           // Debug.WriteLine(IsPointInRectangle(p, corners));
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

            Point3D v1 = new Point3D(corners[1].X - corners[0].X, corners[1].Y - corners[0].Y, corners[1].Z - corners[0].Z);
            Point3D v2 = new Point3D(corners[3].X - corners[0].X, corners[3].Y - corners[0].Y, corners[3].Z - corners[0].Z);

            for (float i = 0; i <= 1; i += 0.005f)
            {
                for (float j = 0; j <= 1; j += 0.005f)
                {
                    Point3D p = new Point3D(
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

        private Point3D ProjectPointToEdge(Point3D p, Point3D[] corners)
        {
            Point3D closestPoint = corners[0];
            double minDist = double.MaxValue;

            for (int i = 0; i < 4; i++)
            {
                Point3D p1 = corners[i];
                Point3D p2 = corners[(i + 1) % 4];
                Point3D projection = ProjectPointOntoLineSegment(p, p1, p2);
                double dist = DistanceSquared(p, projection);

                if (dist < minDist)
                {
                    closestPoint = projection;
                    minDist = dist;
                }
            }

            return closestPoint;
        }

        private Point3D ProjectPointOntoLineSegment(Point3D p, Point3D p1, Point3D p2)
        {
            Point3D v = new Point3D(p2.X - p1.X, p2.Y - p1.Y, p2.Z - p1.Z);
            Point3D w = new Point3D(p.X - p1.X, p.Y - p1.Y, p.Z - p1.Z);

            double c1 = DotProduct(w, v);
            if (c1 <= 0) return p1;

            double c2 = DotProduct(v, v);
            if (c2 <= c1) return p2;

            double b = c1 / c2;
            return new Point3D(
                p1.X + b * v.X,
                p1.Y + b * v.Y,
                p1.Z + b * v.Z
            );
        }

        private double DistanceSquared(Point3D p1, Point3D p2)
        {
            return (p1.X - p2.X) * (p1.X - p2.X) + (p1.Y - p2.Y) * (p1.Y - p2.Y) + (p1.Z - p2.Z) * (p1.Z - p2.Z);
        }


        private double DotProduct(Point3D a, Point3D b)
        {
            return a.X * b.X + a.Y * b.Y + a.Z * b.Z;
        }

        private bool IsPointInRectangle(Point3D p, Point3D[] corners)
        {
            PointF projectedP = new PointF((float)p.X, (float)p.Y);
            PointF[] projectedCorners = new PointF[]
            {
            new PointF((float) corners[0].X, (float) corners[0].Y),
            new PointF((float) corners[1].X, (float) corners[1].Y),
            new PointF((float) corners[2].X, (float) corners[2].Y),
            new PointF((float) corners[3].X, (float) corners[3].Y)
            };

            // Calculate vectors from corners
            Point3D v0 = new Point3D(projectedCorners[3].X - projectedCorners[0].X, projectedCorners[3].Y - projectedCorners[0].Y, 0);
            Point3D v1 = new Point3D(projectedCorners[1].X - projectedCorners[0].X, projectedCorners[1].Y - projectedCorners[0].Y, 0);
            Point3D v2 = new Point3D(projectedP.X - projectedCorners[0].X, projectedP.Y - projectedCorners[0].Y, 0);

            // Calculate dot products
            double dot00 = DotProduct(v0, v0);
            double dot01 = DotProduct(v0, v1);
            double dot02 = DotProduct(v0, v2);
            double dot11 = DotProduct(v1, v1);
            double dot12 = DotProduct(v1, v2);

            // Compute barycentric coordinates
            double invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
            double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
            double v = (dot00 * dot12 - dot01 * dot02) * invDenom;

            // Check if point is in rectangle
            //Debug.WriteLine(u + ' ' + v);
            return (u >= 0) && (v >= 0) && (u <= 1) && (v <= 1);
        }

        private Point3D[] CalculatePlaneCorners(Point3D center, Point3D normal, float width, float height)
        {
            // Calculate two vectors in the plane
            Point3D v1, v2;
            if (normal.X != 0 || normal.Y != 0)
            {
                v1 = new Point3D(-normal.Y, normal.X, 0);
            }
            else
            {
                v1 = new Point3D(0,-normal.Z, normal.Y);
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
            new Point3D(center.X + v1.X + v2.X, center.Y + v1.Y + v2.Y, center.Z + v1.Z + v2.Z),
            new Point3D(center.X + v1.X - v2.X, center.Y + v1.Y - v2.Y, center.Z + v1.Z - v2.Z),
            new Point3D(center.X - v1.X - v2.X, center.Y - v1.Y - v2.Y, center.Z - v1.Z - v2.Z),
            new Point3D(center.X - v1.X + v2.X, center.Y - v1.Y + v2.Y, center.Z - v1.Z + v2.Z)
            };
        }

        private Point3D CrossProduct(Point3D a, Point3D b)
        {
            return new Point3D(
                a.Y * b.Z - a.Z * b.Y,
                a.Z * b.X - a.X * b.Z,
                a.X * b.Y - a.Y * b.X
            );
        }

        private Point3D Normalize(Point3D v)
        {
            float length = (float)Math.Sqrt(v.X * v.X + v.Y * v.Y + v.Z * v.Z);
            return new Point3D(v.X / length, v.Y / length, v.Z / length);
        }

        private Point3D Scale(Point3D v, float scale)
        {
            return new Point3D(v.X * scale, v.Y * scale, v.Z * scale);
        }

        private void DrawAxis(Func<Point3D, Point3D> project, Func<Point3D, Point3D> rotate, Color color, Point3D start, Point3D end)
        {
            using (Pen pen = new Pen(color, 2))
            {
                Point3D rotatedStart = project(rotate(start));
                Point3D rotatedEnd = project(rotate(end));


                graphicsContext.DrawLine(pen, new((int) rotatedStart.X, (int)rotatedStart.Y), new((int)rotatedEnd.X, (int)rotatedEnd.Y));

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
