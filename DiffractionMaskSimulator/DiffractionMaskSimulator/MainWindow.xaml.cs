using System.Windows;
using System.Windows.Controls;

namespace DiffractionMaskSimulator
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        private readonly MainWindowViewModel VM = new MainWindowViewModel();

        public MainWindow()
        {
            InitializeComponent();

            // Set the data context of the window to the view model
            DataContext = VM;
        }

    }
}
