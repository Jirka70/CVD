namespace CVD.launch
{
    internal static class Program
    {

        [STAThread]
        static void Main()
        {
            ApplicationConfiguration.Initialize();
            Application.Run(new VoronoiForm());
        }
    }
}