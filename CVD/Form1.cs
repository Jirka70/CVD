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
        private static readonly int POINT_CLOUD_SIZE = 42;
        private static readonly int NUMBER_OF_ITERATIONS = 10;

        private readonly TextBox xAngleTextBox = new();
        private readonly TextBox yAngleTextBox = new();
        private readonly TextBox zAngleTextBox = new();
        private readonly TextBox scaleBox = new();
        private readonly Button applyButton = new();
        private readonly CheckBox checkBox = new();
        private readonly Button exportBtn = new();

        private readonly Graphics graphicsContext;
        private readonly PictureBox canvasBox = new();
        private readonly List<VoronoiCell3D> diagram;


        public VoronoiForm()
        {
            InitializeComponent();
            canvasBox.Image = new Bitmap(10000, 10000);
            canvasBox.Size = new Size(10000, 10000);
            canvasBox.Location = new Point(0, 0);
            Controls.Add(canvasBox);


            PictureBox exportBox = new();
            exportBox.Size = new Size(100, 100);
            exportBox.Image = new Bitmap(100, 100);

            exportBox.Location = new(WIDTH - 100, 0);
            exportBox.Controls.Add(exportBtn);
            exportBtn.Text = "Export";
            Controls.Add(exportBox);


            Bitmap image = (Bitmap)canvasBox.Image;
            graphicsContext = Graphics.FromImage(image);

            Label xAngleLabel = new();
            Label yAngleLabel = new();
            Label zAngleLabel = new();
            Label textBoxLabel = new();
            Label scaleLabel = new();

            SetUpLabel(xAngleLabel, "X Angle", 0);
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
            Point3D basePoint = new(0, 0, 0);
            Point3D normalVector = new(0, 0, 1);
            Plane2D mainPlane = new(normalVector, basePoint);
            Point3D[] corners = CalculatePlaneCorners(basePoint, normalVector, WIDTH, HEIGHT);


            // Set initial angles and hook up the event
            xAngleTextBox.Text = "0";
            yAngleTextBox.Text = "0";
            zAngleTextBox.Text = "0";
            scaleBox.Text = "1";

            List<Point3D> randomPointCloud = Point3DGenerator.GenerateRandomPoints(POINT_CLOUD_SIZE, WIDTH, HEIGHT);


            diagram = CalculateVoronoiDiagram(randomPointCloud, mainPlane, normalVector, basePoint);
            List<VoronoiCell3D> projectedDiagram = ApplyRotation(diagram, 0, 0, 0, 1, mainPlane, normalVector, basePoint);
            foreach (Point3D cell in randomPointCloud)
            {
                Debug.WriteLine(cell.ToString());
            }
            DrawVoronoiCells(projectedDiagram);
            applyButton.Click += (o, e) => ApplyButton_Click(o, e, mainPlane, normalVector, basePoint, diagram);
        }

        private void Export()
        {
            using (OpenFileDialog openFileDialog = new OpenFileDialog())
            {
                openFileDialog.Filter = "Text files (*.txt)|*.txt|All files (*.*)|*.*";
                openFileDialog.Title = "Vyberte soubor";

                if (openFileDialog.ShowDialog() == DialogResult.OK)
                {
                    //txtFilePath.Text = openFileDialog.FileName;
                }
            }
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

        private void SetUpCheckBox(CheckBox checkBox, int rowNumber)
        {
            int rowHeight = 30;
            checkBox.Location = new(10, rowHeight * rowNumber + 10);
            checkBox.AutoSize = true;
            checkBox.BackColor = Color.Transparent;
            canvasBox.Controls.Add(checkBox);
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
                Vector3 point = new((float)p.X, (float)p.Y, (float)p.Z);
                float radiansX = MathF.PI * angleX / 180f;
                float radiansY = MathF.PI * angleY / 180f;
                float radiansZ = MathF.PI * angleZ / 180f;

                Matrix4x4 rotationX = Matrix4x4.CreateRotationX(radiansX);
                Matrix4x4 rotationY = Matrix4x4.CreateRotationY(radiansY);
                Matrix4x4 rotationZ = Matrix4x4.CreateRotationZ(radiansZ);

                Matrix4x4 combinedRotation = rotationZ * rotationY * rotationX;

                Matrix4x4 scaleMatrix = Matrix4x4.CreateScale(scale);

                Vector3 scaledPoint = Vector3.Transform(point, scaleMatrix);
                Vector3 value = Vector3.Transform(scaledPoint, combinedRotation);
                return new(value.X, value.Y, value.Z);
            };

            Func<Point3D, Point3D> project = (p) =>
            {
                double x = WIDTH / 2 + p.X * scale;
                double y = HEIGHT / 2 - p.Y * scale;
                return new(x, y, 0);
            };

            //DrawPlane(normalVector, basePoint, rotate, project);
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
            //diagram = CenterVoronoiDiagram(diagram, NUMBER_OF_ITERATIONS);
            List<VoronoiCell3D> projectedDiagram = ProjectVoronoiCellsReversively(diagram, pointTransformation);

            return projectedDiagram;
        }

        private List<VoronoiCell3D> ProjectVoronoiCellsOnPlane(Plane2D plane, Point3D basePoint, Func<Point3D, Point3D> rotate, Func<Point3D, Point3D> project, Point3D normalVector, List<VoronoiCell3D> voronoiCells)
        {
            Point3D[] corners = CalculatePlaneCorners(basePoint, normalVector, WIDTH, HEIGHT);
            Debug.WriteLine(corners[0].ToString() + " " + corners[1].ToString() + " " + corners[2].ToString() + " " + corners[3].ToString());

            Function planeEquation = plane.equation;

            List<VoronoiCell3D> projectedCells = new();

            foreach (VoronoiCell3D voronoiCell in voronoiCells)
            {
                VoronoiCell3D? projectedCell = ProjectVoronoiCell(basePoint, normalVector, corners, rotate, project, voronoiCell);
                if (projectedCell != null)
                {
                    projectedCells.Add(projectedCell);
                }
            }

            DrawAxis(project, rotate, Color.Red, new Point3D(-500, 0, 0), new Point3D(500, 0, 0), "X");
            DrawAxis(project, rotate, Color.Green, new Point3D(0, -500, 0), new Point3D(0, 500, 0), "Y");
            DrawAxis(project, rotate, Color.Blue, new Point3D(0, 0, -500), new Point3D(0, 0, 500), "Z");
            //DrawPlane(normalVector, basePoint, rotate, project);

            return projectedCells;
        }

        private VoronoiCell3D? ProjectVoronoiCell(Point3D basePoint, Point3D normalVector, Point3D[] corners, Func<Point3D, Point3D> rotate, Func<Point3D, Point3D> project, VoronoiCell3D voronoiCell)
        {
            Point3D center = voronoiCell.getCenter();
            List<Edge3D> edges = voronoiCell.GetEdges();

            List<Edge3D> projectedEdges = new();
            foreach (Edge3D edge in edges)
            {
                ProjectedEdge? projectedEdge = ProjectVoronoiEdge(basePoint, normalVector, corners, rotate, project, edge);
                if (projectedEdge == null)
                {
                    continue;
                }

                if (projectedEdge.InScreenEdge != null)
                {
                    projectedEdges.Add(projectedEdge.InScreenEdge);
                }

                if (projectedEdge.OutScreenEdge != null)
                {
                    projectedEdges.Add(projectedEdge.OutScreenEdge);
                }
            }

            Point3D projectedCenter = RotatePoint(basePoint, normalVector, corners, rotate, project, center);

            return new(project(projectedCenter), projectedEdges);

        }

        private bool IsPointInsideRectangle(Point3D p, Point3D[] corners)
        {
            PointF projectedP = new PointF((float)p.X, (float)p.Y);
            PointF[] projectedCorners = new PointF[]
            {
            new PointF((float) corners[0].X, (float) corners[0].Y),
            new PointF((float) corners[1].X, (float) corners[1].Y),
            new PointF((float) corners[2].X, (float) corners[2].Y),
            new PointF((float) corners[3].X, (float) corners[3].Y)
            };

            Vector2 v0 = new Vector2(projectedCorners[3].X - projectedCorners[0].X, projectedCorners[3].Y - projectedCorners[0].Y);
            Vector2 v1 = new Vector2(projectedCorners[1].X - projectedCorners[0].X, projectedCorners[1].Y - projectedCorners[0].Y);
            Vector2 v2 = new Vector2(projectedP.X - projectedCorners[0].X, projectedP.Y - projectedCorners[0].Y);

            float dot00 = DotProduct(v0, v0);
            float dot01 = DotProduct(v0, v1);
            float dot02 = DotProduct(v0, v2);
            float dot11 = DotProduct(v1, v1);
            float dot12 = DotProduct(v1, v2);

            float invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
            float u = (dot11 * dot02 - dot01 * dot12) * invDenom;
            float v = (dot00 * dot12 - dot01 * dot02) * invDenom;

            return (u >= 0) && (v >= 0) && (u <= 1) && (v <= 1);
        }

        private float DotProduct(Vector2 a, Vector2 b)
        {
            return a.X * b.X + a.Y * b.Y;
        }

        private Point3D RotatePoint(Point3D basePoint, Point3D normalVector, Point3D[] corners, Func<Point3D, Point3D> rotate, Func<Point3D, Point3D> project, Point3D pointToProject)
        {
            double z = (normalVector.X * (pointToProject.X - basePoint.X) + normalVector.Y * (pointToProject.Y - basePoint.Y)) / -normalVector.Z + basePoint.Z;
            Point3D p = new(pointToProject.X, pointToProject.Y, z);
            Point3D rotatedPoint = rotate(p);

            return rotatedPoint;
        }

        private ProjectedEdge? ProjectVoronoiEdge(Point3D basePoint, Point3D normalVector, Point3D[] corners, Func<Point3D, Point3D> rotate, Func<Point3D, Point3D> project, Edge3D edgeToProject)
        {
            Point3D start = edgeToProject.startingPoint;
            Point3D end = edgeToProject.endingPoint;

            Point3D rotatedStart = RotatePoint(basePoint, normalVector, corners, rotate, project, start);
            Point3D rotatedEnd = RotatePoint(basePoint, normalVector, corners, rotate, project, end);

            bool isStartInRectangle = IsPointInsideRectangle(rotatedStart, corners);
            bool isEndInRectangle = IsPointInsideRectangle(rotatedEnd, corners);

            if (isStartInRectangle && isEndInRectangle)
            {
                return new(new(project(rotatedStart), project(rotatedEnd)), null);

            }

            if (!isStartInRectangle && !isEndInRectangle)
            {
                return new(null, new(project(rotatedStart), project(rotatedEnd)));
            }

            if (!isStartInRectangle) // end is contained within the rectangle
            {
                Point3D unprojectedRotatedStart = rotatedStart;
                return new(new(project(rotatedEnd), project(rotatedEnd)),
                    new(project(rotatedEnd), project(unprojectedRotatedStart)));
            }

            if (!isEndInRectangle) // start is contained within the rectangle
            {
                Point3D unprojectedRotatedEnd = rotatedEnd;
                return new(new(project(rotatedStart), project(rotatedEnd)), new(project(rotatedEnd), project(unprojectedRotatedEnd)));
            }

            return null;
        }

        private Point3D FindIntersection(Point3D lineStart, Point3D lineEnd, Point3D edgeStart, Point3D edgeEnd)
        {
            Vector2 lineDir = new Vector2((float)(lineEnd.X - lineStart.X), ((float)(lineEnd.Y - lineStart.Y)));
            Vector2 edgeDir = new Vector2((float)(edgeEnd.X - edgeStart.X), (float)(edgeEnd.Y - edgeStart.Y));
            Vector2 lineToEdgeStart = new Vector2((float)(edgeStart.X - lineStart.X), ((float)(edgeStart.Y - lineStart.Y)));

            float a = DotProduct(lineDir, lineDir);
            float b = DotProduct(lineDir, edgeDir);
            float e = DotProduct(edgeDir, edgeDir);
            float d = a * e - b * b;

            if (Math.Abs(d) < 1e-6)
            {
                // Lines are parallel
                return null;
            }

            float c = DotProduct(lineDir, lineToEdgeStart);
            float f = DotProduct(edgeDir, lineToEdgeStart);

            float t = (b * f - c * e) / d;
            float u = (a * f - b * c) / d;

            if (t < 0 || t > 1 || u < 0 || u > 1)
            {
                // Intersection point is not within the line segment or the edge
                return null;
            }

            return new Point3D(
                lineStart.X + t * lineDir.X,
                lineStart.Y + t * lineDir.Y,
                0
            );
        }

        private Point3D? GetIntersectionPoint(Point3D lineStart, Point3D lineEnd, Point3D[] corners)
        {
            Point3D? closestIntersection = null;
            double minDistance = double.MaxValue;

            var edges = new[]
            {
        new { Start = corners[0], End = corners[1] },
        new { Start = corners[1], End = corners[2] },
        new { Start = corners[2], End = corners[3] },
        new { Start = corners[3], End = corners[0] }
    };

            foreach (var edge in edges)
            {
                Point3D intersection = FindIntersection(lineStart, lineEnd, edge.Start, edge.End);

                if (intersection != null)
                {
                    double distance = GetDistance(lineStart, intersection);

                    if (distance < minDistance)
                    {
                        minDistance = distance;
                        closestIntersection = intersection;
                    }
                }
            }

            return closestIntersection;
        }

        private double GetDistance(Point3D point1, Point3D point2)
        {
            return Math.Sqrt(
                Math.Pow(point1.X - point2.X, 2) +
                Math.Pow(point1.Y - point2.Y, 2) +
                Math.Pow(point1.Z - point2.Z, 2)
            );
        }


        private double DistanceSquared(Point3D p1, Point3D p2)
        {
            return (p1.X - p2.X) * (p1.X - p2.X) + (p1.Y - p2.Y) * (p1.Y - p2.Y) + (p1.Z - p2.Z) * (p1.Z - p2.Z);
        }

        private void DrawPlane(Point3D normal, Point3D point, Func<Point3D, Point3D> rotate, Func<Point3D, Point3D> project)
        {
            Point3D[] corners = CalculatePlaneCorners(point, normal, WIDTH, HEIGHT);
            Plane2D plane = new(normal, point);

            Debug.WriteLine(corners[0].ToString() + " " + corners[1].ToString() + " " + corners[2].ToString() + " " + corners[3].ToString());
            Debug.WriteLine(IsPointInsideRectangle(new(400, 200, 0), corners));

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
                        graphicsContext.FillRectangle(semiTransparentBrush, (int)p2D.X, (int)p2D.Y, 2, 2); // Draw point with semi-transparent color
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
                v1 = new(0, -normal.Z, normal.Y);
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
            new(center.X + v1.X + v2.X, center.Y + v1.Y + v2.Y, 0),
            new(center.X + v1.X - v2.X, center.Y + v1.Y - v2.Y, 0),
            new(center.X - v1.X - v2.X, center.Y - v1.Y - v2.Y, 0),
            new(center.X - v1.X + v2.X, center.Y - v1.Y + v2.Y, 0)
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
            graphicsContext.DrawLine(new Pen(color, 2), s, e);

            DrawAxisLabel(project, rotate, color, start, axisName);
            DrawAxisLabel(project, rotate, color, end, axisName);

            int interval = 100;
            for (int i = (int)(start.X + interval); i <= end.X - interval; i += interval)
            {
                if (axisName == "X")
                    DrawAxisLabel(project, rotate, color, new Point3D(i, start.Y, start.Z), i.ToString());
            }
            for (int i = (int)(start.Y + interval); i <= end.Y - interval; i += interval)
            {
                if (axisName == "Y")
                    DrawAxisLabel(project, rotate, color, new Point3D(start.X, i, start.Z), i.ToString());
            }
            for (int i = (int)(start.Z + interval); i <= end.Z - interval; i += interval)
            {
                if (axisName == "Z")
                    DrawAxisLabel(project, rotate, color, new Point3D(start.X, start.Y, i), i.ToString());
            }
        }

        private void DrawAxisLabel(Func<Point3D, Point3D> project, Func<Point3D, Point3D> rotate, Color color, Point3D position, string labelValue)
        {
            Point3D projectedPosition = project(rotate(position));
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
            Point3D groundNormalVector = new(0, 0, 1);
            Point3D initialPoint = new(0, 0, 0);
            Point3D normalVector = plane.normalVector;
            Plane2D groundPlane = new(groundNormalVector, initialPoint);

            Point3D rotationAxisVector = normalVector.CalculateCrossProduct(groundNormalVector);
            double angle = plane.CalculateAngle(groundPlane);

            Rotation rotation = new(rotationAxisVector, (float)angle);
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
                }
                else
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
            double y = (-d - a * x) / b;
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
                    new((int)vertex.startingPoint.X, (int)vertex.startingPoint.Y),
                    new((int)vertex.endingPoint.X, (int)vertex.endingPoint.Y));
            }

            Draw3DPoint(pointPen, cell.getCenter());
        }

        private class ProjectedEdge
        {
            public readonly Edge3D? InScreenEdge;
            public readonly Edge3D? OutScreenEdge;

            public ProjectedEdge(Edge3D? inScreenEdge, Edge3D? outScreenEdge)
            {
                InScreenEdge = inScreenEdge;
                OutScreenEdge = outScreenEdge;
            }
        }

        private class Line
        {
            public readonly Point3D p1;
            public readonly Point3D p2;

            public readonly Vector3 equation;

            public Line(Point3D point1, Point3D point2)
            {
                p1 = point1;
                p2 = point2;
                Vector2 vec = new((float)Math.Abs(p1.X - p2.X), (float)Math.Abs(p1.Y - p2.Y));
                Vector2 normal = new(-vec.X, vec.Y);

                float c = (float)(-normal.X * p1.X - normal.Y * p1.Y);
                equation = new(normal.X, normal.Y, c);
            }
        }

        private void button1_Click(object sender, EventArgs e)
        {

        }
    }
}
