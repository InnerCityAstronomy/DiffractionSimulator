using System;
using System.Drawing;
using System.Linq;
using System.Numerics;
using System.Runtime.InteropServices;
using System.Threading.Tasks;

namespace DiffractionMaskSimulator
{
    public class MainWindowModel
    {
        // Public variables
        /// <summary> Filename of the mask image.</summary>
        public string Filename
        {
            get
            {
                return _Filename;
            }
            set
            {
                _Filename = value;
                MaskLoaded = false;
            }
        }
        
        /// <summary> Diameter of the telescope aperture, in metres.</summary>
        public double DAperture { get; set; }
        
        /// <summary> Focal length of the telescope, in metres.</summary>
        public double FL { get; set; }

        /// <summary> Smallest wavelength in the wavelengths array in metres.</summary>
        public double LambdaStart { get; set; }

        /// <summary> Number of wavelengths in the wavelengths array.</summary>
        public int NumLambda { get; set; }

        /// <summary> Largest wavelength in the wavelengths array in metres.</summary>
        public double LambdaStop { get; set; }
        
        /// <summary> The focus error of the telescope, in metres. A value of zero corresponds to perfect focus.</summary>
        public double Defocus { get; set; }

        /// <summary> Brightness added to the final image, in decibels.</summary>
        public double BrightnessdB { get; set; }

        /// <summary> Resolution of the diffraction image, in pixels.</summary>
        public int Resolution
        {
            get
            {
                return _Resolution;
            }
            set
            {
                _Resolution = value;
                MaskLoaded = false;
            }
        }
        
        /// <summary> If set to true, the amplitudes of the wavelengths
        /// are scaled according to the D65 standard illuminant .</summary>
        public bool UseD65Weighting { get; set; }

        // Private Variables
        private string _Filename { get; set; }
        private int _Resolution { get; set; }
        private double LMaskMirror { get; set; }
        private bool MaskLoaded { get; set; }
        private Image MaskImage { get; set; }
        private Complex[] DiffractionMaster { get; set; }
        private Bitmap DiffractionPattern { get; set; }

        /// <summary> If set to true, will use all available cores for the computation.</summary>
        private bool UseMultiThreading { get; set; }
        private double dx { get; set; }
        private int[] Dimensions { get; set; }

        // Settings
        private const double ImageDynamicRange = 24;
        private const double Gamma = 1.0;
        private MathNet.Numerics.IntegralTransforms.FourierOptions FOptions { get; set; }

        // Default Constructor
        public MainWindowModel() 
        {
            // Public parameters
            Filename = "Bahtinov.jpg";
            DAperture = 0.2;
            FL = 1.2;
            Defocus = 0E-3;
            LambdaStart = 400E-9;
            LambdaStop = 700E-9;
            NumLambda = 13;
            Resolution = 1200;
            BrightnessdB = 30;
            UseD65Weighting = true;

            // Private parameters
            LMaskMirror = 1;            
            MaskLoaded = false;
            UseMultiThreading = true;
            FOptions = PhysicalOptics.FOptions;
        }

        private void LoadMask()
        {
            // Load the mask file
            MaskImage = Image.FromFile(Filename);

            // Resize the image to the desired resolution
            MaskImage = Util.ResizeImage(MaskImage, Resolution, Resolution);
            Dimensions = new int[] { Resolution, Resolution };

            // Convert it to a matrix
            DiffractionMaster = Array.ConvertAll(GetBimtapBrightnessMatrix(new Bitmap(MaskImage)), item => (Complex)item);

            // Compute the master diffraction pattern
            MathNet.Numerics.IntegralTransforms.Fourier.Forward2D(DiffractionMaster, Dimensions[0], Dimensions[1], FOptions);
            DiffractionMaster = Util.FFTShift(DiffractionMaster, Dimensions);

            // Normalise the master diffraction pattern
            double ScaleFac = 1 / Array.ConvertAll(DiffractionMaster, item => Complex.Abs(item)).Max();

            for (int i = 0; i < DiffractionMaster.Length; i++)
            {
                DiffractionMaster[i] *= ScaleFac;
            }

            // Set ImageLoaded to true so we don't have to load it every time
            MaskLoaded = true;
        }

        /// <summary>
        /// Will check all of the settings and throw an <see cref="Exception"/> if any are invalid.
        /// </summary>
        /// <returns><see cref="void"/></returns>
        private void ValidateSettings()
        {
            if (Resolution <= 1)
            {
                throw new ArgumentException(string.Format("The resolution must be greater than 1."));
            }
            else if (DAperture <= 0)
            {
                throw new ArgumentException(string.Format("The aperture diameter must be greater than zero."));
            }
            else if (FL <= 0)
            {
                throw new ArgumentException(string.Format("The focal length must be greater than zero."));
            }
            else if (NumLambda <= 0)
            {
                throw new ArgumentException(string.Format("The number of wavelength points must be greater than zero."));
            }
            else if (LambdaStop <= LambdaStart)
            {
                // Swap the start and stop wavelengths
                double TEMP = LambdaStart;
                LambdaStart = LambdaStop;
                LambdaStop = TEMP;
            }
        }

        /// <summary>
        /// Main function which computes the diffraction pattern.
        /// </summary>
        /// <returns>The diffraction pattern as an <see cref="Image"/></returns>
        public Image GenerateImage()
        {
            // Check that the settings are valid
            ValidateSettings();

            // Construct the wavelengths array
            double[] Wavelengths = MathNet.Numerics.Generate.LinearSpaced(NumLambda, LambdaStart, LambdaStop);

            // If the mask hasn't been loaded, load it
            if (!MaskLoaded)
            {
                LoadMask();
            }

            // Compute the x and y meshgrid as this only needs to be calculated once
            double[] xvals = MathNet.Numerics.Generate.LinearSpaced(Resolution, -DAperture / 2.0, DAperture / 2.0);

            // Find the final step sizes
            dx = xvals[1] - xvals[0];

            // Initialise the RGB master data
            int[] RMaster = new int[DiffractionMaster.Length];
            int[] GMaster = new int[DiffractionMaster.Length];
            int[] BMaster = new int[DiffractionMaster.Length];

            // Compute the smallest wavelength first to establish the master XX and YY domains
            (double[] XXMaster, double[] YYMaster, _) = ComputeDiffraction(Wavelengths[0]);

            if (UseMultiThreading)
            {
                Parallel.For(0, Wavelengths.Length, i =>
                {
                    (double[] XX, double[] YY, double[] DiffractionMag) = ComputeDiffraction(Wavelengths[i]);

                    if (i != 0)
                    {
                        double[] x = XX.Distinct().ToArray();
                        double[] y = YY.Distinct().ToArray();

                        // Interpolate the magnitude data to the master grid
                        DiffractionMag = Util.BicubicInterpolation(x, y, Util.Reshape1DTo2D(DiffractionMag, Dimensions), XXMaster, YYMaster);
                    }

                    // Apply D65 amplitude weighting
                    if (UseD65Weighting)
                    {
                        DiffractionMag = ApplyD65Weighting(DiffractionMag, Wavelengths[i]);
                    }

                    // Scale the data for the image
                    DiffractionMag = ScaleData(DiffractionMag);

                    Color Col = PhysicalOptics.WavelengthToColour(Wavelengths[i], 1, Gamma);

                    lock (RMaster)
                    {
                        lock (GMaster)
                        {
                            lock (BMaster)
                            {
                                // Add the colourdata to the RGB master components
                                for (int j = 0; j < DiffractionMag.Length; j++)
                                {
                                    double TEMP_MAG = DiffractionMag[j];
                                    double TEMP_R = Col.R;
                                    double TEMP_RMASTER = (int)(TEMP_R * TEMP_MAG);

                                    RMaster[j] += (int)(Col.R * DiffractionMag[j]);
                                    GMaster[j] += (int)(Col.G * DiffractionMag[j]);
                                    BMaster[j] += (int)(Col.B * DiffractionMag[j]);
                                }
                            }
                        }
                    }

                });
            }
            else
            {
                for (int i = 0; i < Wavelengths.Length; i++)
                {
                    (double[] XX, double[] YY, double[] DiffractionMag) = ComputeDiffraction(Wavelengths[i]);

                    if (i != 0)
                    {
                        double[] x = XX.Distinct().ToArray();
                        double[] y = YY.Distinct().ToArray();

                        // Interpolate the magnitude data to the master grid
                        DiffractionMag = Util.BicubicInterpolation(x, y, Util.Reshape1DTo2D(DiffractionMag, Dimensions), XXMaster, YYMaster);
                    }

                    // Apply D65 amplitude weighting
                    if (UseD65Weighting)
                    {
                        DiffractionMag = ApplyD65Weighting(DiffractionMag, Wavelengths[i]);
                    }

                    // Scale the data for the image
                    DiffractionMag = ScaleData(DiffractionMag);

                    // Add the colourdata to the RGB master components
                    Color Col = PhysicalOptics.WavelengthToColour(Wavelengths[i], 1, Gamma);

                    for (int j = 0; j < DiffractionMag.Length; j++)
                    {
                        RMaster[j] += (int)(Col.R * DiffractionMag[j]);
                        GMaster[j] += (int)(Col.G * DiffractionMag[j]);
                        BMaster[j] += (int)(Col.B * DiffractionMag[j]);
                    }
                }
            }

            // Ensure RGB data is between 0 and 255
            for (int i = 0; i < RMaster.Length; i++)
            {
                RMaster[i] = RMaster[i] / Wavelengths.Length;
                GMaster[i] = GMaster[i] / Wavelengths.Length;
                BMaster[i] = BMaster[i] / Wavelengths.Length;

                RMaster[i] = RMaster[i].Clamp(0, 255);
                GMaster[i] = GMaster[i].Clamp(0, 255);
                BMaster[i] = BMaster[i].Clamp(0, 255);
            }            

            // Combine the RGB components
            Color[] ColourMaster = CombineRGBComponents(RMaster, GMaster, BMaster);

            // Convert to a bitmap
            DiffractionPattern = ColourDataToBitmap(ColourMaster, Dimensions[1], Dimensions[0]);

            //return MaskImage;
            return DiffractionPattern;
        }

        /// <summary>
        /// Computes the magnitude of the diffraction pattern at the prime focus of the telescope.
        /// </summary>
        /// <param name="lambda">The wavelength of the light.</param>
        /// <returns>A tuple containing:<para/>
        /// The X values in meshgrid format<para/>
        /// The Y values in meshgrid format<para/>
        /// The magnitude of the diffraction pattern</returns>
        public (double[], double[], double[]) ComputeDiffraction(double lambda)
        {
            (double[] x, double[] y) = GetDiffractionDomain(lambda);
            (double[] XX, double[] YY) = Util.MeshGrid1D(x, y);

            // Clone the master diffraction pattern
            Complex[] Diffraction = (Complex[])DiffractionMaster.Clone();

            // Apply the lens/primary mirror transformation to the field
            Diffraction = PhysicalOptics.ApplyThinLensArray(XX, YY, Diffraction, lambda, FL);

            // Propagate the field by the defocus amount
            if (Defocus != 0)
            {
                (XX, YY, Diffraction) = PhysicalOptics.PropagateField(XX, YY, Diffraction, Dimensions, lambda, Defocus);
            }

            // Extract the intensity of the diffraction pattern, which is the square of the E field
            double[] DiffractionMag = Array.ConvertAll(Diffraction, item => Complex.Abs(item) * Complex.Abs(item));

            return (XX, YY, DiffractionMag);
        }

        /// <summary>
        /// Computes the x and y values of the transformed region.
        /// </summary>
        /// <param name="lambda">The wavelength of the light.</param>
        /// <returns>A tuple containing the unique x and y values.</returns>
        public (double[], double[]) GetDiffractionDomain(double lambda)
        {
            double[] Fx = Util.GetFreqs(Resolution, 1 / dx);
            double[] x = new double[Fx.Length];

            for (int i = 0; i < Fx.Length; i++)
            {
                x[i] = Fx[i] * lambda * LMaskMirror;
            }

            return (x, x);
        }

        /// <summary>
        /// Takes an array of datapoints and scales it between 0 and 255 for use in an image.
        /// </summary>
        /// <param name="data">Array of data arranged in row-major order.</param>
        /// <param name="ImageDynamicRange">The dynamic range (in dB) of the output data.</param>
        /// <returns>The data scaled between 0 and 255.</returns>
        private double[] ScaleData(double[] data)
        {
            // Loop through and scale the data
            for (int i = 0; i < data.Length; i++)
            {
                // The input data is always between 0 and 1
                // Taking 10*Log10 converts to a logarithmic scale and scales data between -Inf and 0
                // Dividing by the dynamic range and adding 1 means that the dynamic range is captured between 0 and 1
                data[i] = (10.0 * Math.Log10(data[i]) + BrightnessdB) / ImageDynamicRange + 1.0;

                // Ignore any result that is NaN
                if (double.IsNaN(data[i]))
                {
                    data[i] = 0;
                }

            }

            return data;
        }

        /// <summary>
        /// Converts a bitmap image to an array of values containing the brightness arranged in row-major order.<para/>
        /// Source: https://stackoverflow.com/questions/13481558/converting-a-bitmap-image-to-a-matrix
        /// </summary>
        /// <param name="b1">The <see cref="Bitmap"/> Image to convert.</param>
        /// <returns>An array containing the brightness of each pixel.</returns>
        private static double[] GetBimtapBrightnessMatrix(Bitmap b1)
        {
            int height = b1.Height;
            int width = b1.Width;
            int[] dims = new[] { b1.Height, b1.Width };

            double[] BrightnessMatrix = new double[width*height];
            for (int i = 0; i < dims[0]; i++)
            {
                for (int j = 0; j < dims[1]; j++)
                {
                    int idx = j + i * dims[1];
                    BrightnessMatrix[idx] = b1.GetPixel(j, i).GetBrightness();
                }
            }
            return BrightnessMatrix;
        }

        /// <summary>
        /// Takes an array of datapoints and applies the D65 standard illuminant amplitude weighting.
        /// </summary>
        /// <param name="data">Array of data arranged in row-major order.</param>
        /// <param name="lambda">Wavelength of light, nm.</param>
        /// <returns>The data scaled between 0 and 255.</returns>
        private static double[] ApplyD65Weighting(double[] data, double lambda)
        {
            double Weight = Util.CubicInterpolation(PhysicalOptics.D65Freqs, PhysicalOptics.D65WeightsNormalised, lambda*1E9);

            for (int i = 0; i < data.Length; i++)
            {
                data[i] *= Weight;
            }

            return data;
        }

        /// <summary>
        /// Converts an array of <see cref="Color"/> objects arranged in row-major order to a <see cref="Bitmap"/> image.
        /// </summary>
        /// <param name="data">Array of <see cref="Color"/> objects arranged in row-major order.</param>
        /// <param name="width">The width of the image.</param>
        /// <param name="height">The height of the image.</param>
        /// <returns>The bitmap image.</returns>
        private static Bitmap ColourDataToBitmap(Color[] data, int width, int height)
        {
            int[] dims = new int[] { height, width };

            Int32[] bits = new Int32[width * height];
            GCHandle handle = GCHandle.Alloc(bits, GCHandleType.Pinned);
            Bitmap bmp = new Bitmap(width, height, width * 4,
              System.Drawing.Imaging.PixelFormat.Format32bppPArgb, handle.AddrOfPinnedObject());

            for (int i = 0; i < dims[0]; i++)
            {
                for (int j = 0; j < dims[1]; j++)
                {
                    int idx = j + i * dims[1];
                    bits[idx] = data[idx].ToArgb();
                }
            }

            return bmp;
        }

        /// <summary>
        /// Combines separate RGB components into a single <see cref="Color"/> object.
        /// </summary>
        /// <param name="RMaster">Array of Red values.</param>
        /// <param name="GMaster">Array of Green values.</param>
        /// <param name="BMaster">Array of Blue values.</param>
        /// <returns>An array of <see cref="Color"/> objects.</returns>
        private static Color[] CombineRGBComponents(int[] RMaster, int[] GMaster, int[] BMaster)
        {
            Color[] Result = new Color[RMaster.Length];

            for (int i = 0; i < Result.Length; i++)
            {
                Result[i] = Color.FromArgb(RMaster[i], GMaster[i], BMaster[i]);
            }

            return Result;
        }

    }
}
