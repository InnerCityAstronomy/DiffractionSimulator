using Microsoft.Win32;
using System.Windows;
using System.Windows.Controls;


namespace DiffractionMaskSimulator
{
    /// <summary>
    /// Interaction logic for SaveFileControl.xaml
    /// </summary>
    public partial class SaveFileControl : UserControl
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
            typeof(SaveFileControl),
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
            typeof(SaveFileControl),
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
            typeof(SaveFileControl),
            new PropertyMetadata()
        );

        public SaveFileControl()
        {
            FilePath = "";

            InitializeComponent();
        }

        private void BrowseButton_Click(object sender, RoutedEventArgs e)
        {
            // Open file dialog window
            SaveFileDialog saveFileDialog = new SaveFileDialog
            {
                Filter = Filter,
                RestoreDirectory = true,
                Title = Title
            };

            if (saveFileDialog.ShowDialog() == true)
            {
                FilePath = saveFileDialog.FileName;
            }

        }
    }
}
