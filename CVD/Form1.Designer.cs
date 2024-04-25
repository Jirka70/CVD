namespace CVD
{
    public partial class VoronoiForm : Form
    {
        private static readonly int WIDTH = 800;
        private static readonly int HEIGHT = 600;

        private System.ComponentModel.IContainer components = null;


        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows Form Designer generated code

        /// <summary>
        ///  Required method for Designer support - do not modify
        ///  the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            this.components = new System.ComponentModel.Container();
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(WIDTH, HEIGHT);
            this.Text = "Voronoi Diagrams - Jiri Tresohlavy";
        }

        #endregion
    }
}
