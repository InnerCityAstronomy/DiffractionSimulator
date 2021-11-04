using System.Drawing;
using System.Windows.Input;
using System.ComponentModel;
using Microsoft.Win32;
using System;
using System.Windows;

namespace DiffractionMaskSimulator
{
    public class MainWindowViewModel : INotifyPropertyChanged
    {
        // Model
        private MainWindowModel Model = new MainWindowModel();

        // Variables
        public string Filename
        {
            get
            {
                return Model.Filename;
            }
            set
            {
                Model.Filename = value;
                OnPropertyChange("Filename");
            }
        }
        public double DAperturemm
        {
            get
            {
                return Model.DAperture*1E3;
            }
            set
            {
                Model.DAperture = value * 1E-3;
                OnPropertyChange("DAperturemm");
            }
        }
        public double FLmm
        {
            get
            {
                return Model.FL * 1E3;
            }
            set
            {
                Model.FL = value * 1E-3;
                OnPropertyChange("FLmm");
            }
        }
        public double Defocusmm
        {
            get
            {
                return Model.Defocus * 1E3;
            }
            set
            {
                Model.Defocus = value * 1E-3;
                OnPropertyChange("Defocusmm");
            }
        }
        public double LambdaStartnm
        {
            get
            {
                return Model.LambdaStart * 1E9;
            }
            set
            {
                Model.LambdaStart = value * 1E-9;
                OnPropertyChange("LambdaStartnm");
            }
        }
        public double LambdaStopnm
        {
            get
            {
                return Model.LambdaStop * 1E9;
            }
            set
            {
                Model.LambdaStop = value * 1E-9;
                OnPropertyChange("LambdaStopnm");
            }
        }
        public int NumLambda
        {
            get
            {
                return Model.NumLambda;
            }
            set
            {
                Model.NumLambda = value;
                OnPropertyChange("NumLambda");
            }
        }
        public double BrightnessdB
        {
            get
            {
                return Model.BrightnessdB;
            }
            set
            {
                Model.BrightnessdB = value;
                OnPropertyChange("BrightnessdB");
            }
        }
        public int Resolution
        {
            get
            {
                return Model.Resolution;
            }
            set
            {
                Model.Resolution = value;
                OnPropertyChange("Resolution");
            }
        }
        public bool UseD65Weighting
        {
            get
            {
                return Model.UseD65Weighting;
            }
            set
            {
                Model.UseD65Weighting = value;
                OnPropertyChange("UseD65Weighting");
            }
        }
        public bool CanSaveImage
        {
            get
            {
                return ImageMask != null;
            }
        }
        public bool CanExecute
        {
            get
            {
                // check if executing is allowed, i.e., validate, check if a process is running, etc. 
                return true;
            }
        }

        // Image
        private Image _ImageMask { get; set; }
        public Image ImageMask
        {
            get
            {
                return _ImageMask;
            }
            set
            {
                _ImageMask = value;
                OnPropertyChange("ImageMask");
                OnPropertyChange("CanSaveImage");
            }
        }

        // Commands
        private ICommand _CmdGenImage;
        public ICommand CmdGenImage
        {
            get
            {
                return _CmdGenImage ??= new CommandHandler(() =>
                {
                    try
                    {
                        ImageMask = Model.GenerateImage();
                    }
                    catch (Exception ex)
                    {
                        MessageBox.Show(ex.Message);
                    }
                    
                }, () => CanExecute);

            }
        }

        private ICommand _CmdSaveImage;
        public ICommand CmdSaveImage
        {
            get
            {
                return _CmdSaveImage ??= new CommandHandler(() =>
                {
                    SaveFileDialog SFDialog = new SaveFileDialog();
                    SFDialog.Filter = "JPEG|*.jpg|Bitmap|*.bmp|GIF|*.gif|PNG|*.png";
                    SFDialog.Title = "Save image";

                    if (SFDialog.ShowDialog() == true)
                    {
                        ImageMask.Save(SFDialog.FileName);
                    }
                    
                }, () => CanExecute);
            }
        }

        public MainWindowViewModel()
        {
        }

        // Implement required methods for INotifyPropertyChanged
        public event PropertyChangedEventHandler PropertyChanged;

        protected void OnPropertyChange(string propertyName)
        {
            if (PropertyChanged != null)
            {
                PropertyChanged(this, new PropertyChangedEventArgs(propertyName));
            }
        }

    }
}
