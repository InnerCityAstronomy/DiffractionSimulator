using Microsoft.Win32;
using System.Windows;
using System.Windows.Controls;


namespace DiffractionMaskSimulator
{
    /// <summary>
    /// Interaction logic for BrowseFileControl.xaml
    /// </summary>
    public partial class BrowseFileControl : UserControl
    {
        public string Title
        {
            get { return (string)GetValue(TitleProperty); }
            set { SetValue(TitleProperty, value); }
        }
        public static readonly DependencyProperty TitleProperty =
        DependencyProperty.Register
        (
            "Title",
            typeof(string),
            typeof(BrowseFileControl),
            new PropertyMetadata()
        );

        public string Filter
        {
            get { return (string)GetValue(FilterProperty); }
            set { SetValue(FilterProperty, value); }
        }
        public static readonly DependencyProperty FilterProperty =
        DependencyProperty.Register
        (
            "Filter",
            typeof(string),
            typeof(BrowseFileControl),
            new PropertyMetadata()
        );


        public string FilePath
        {
            get { return (string)GetValue(FilePathProperty); }
            set { SetValue(FilePathProperty, value); }
        }
        public static readonly DependencyProperty FilePathProperty =
        DependencyProperty.Register
        (
            "FilePath",
            typeof(string),
            typeof(BrowseFileControl),
            new PropertyMetadata()
        );

        public BrowseFileControl()
        {
            FilePath = "";

            InitializeComponent();
        }

        private void BrowseButton_Click(object sender, RoutedEventArgs e)
        {
            // Open file dialog window
            OpenFileDialog openFileDialog = new OpenFileDialog
            {
                Filter = Filter,
                RestoreDirectory = true,
                Title = Title
            };

            if (openFileDialog.ShowDialog() == true)
            {
                FilePath = openFileDialog.FileName;
            }

        }
    }
}
