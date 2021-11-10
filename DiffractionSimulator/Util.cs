using System;
using System.Drawing;
using System.Drawing.Drawing2D;
using System.Numerics;

namespace DiffractionMaskSimulator
{
    /// <summary>
    /// A <c>static</c> class which provides a range of useful methods.
    /// </summary>
    public static class Util
    {
        /// <summary> Computes the product of the elements of the <paramref name="Dimensions"/>
        /// array starting from the index specified by <paramref name="StartIdx"/>.</summary>
        /// <param name="Dimensions"> Array which contains the dimensions.</param>
        /// <param name="StartIdx"> The index in <paramref name="Dimensions"/> to start the product from.</param>
        /// <returns> The product.</returns>
        private static int DimensionProduct(int[] Dimensions, int StartIdx)
        {
            int Product = 1;

            for (int i = StartIdx; i < Dimensions.Length; i++)
            {
                Product *= Dimensions[i];
            }

            return Product;
        }

        /// <summary> Converts multidimensional indices to a single equivalent linear index.</summary>
        /// <param name="indices"> The multidimensional indices to convert.</param>
        /// <param name="Dimensions"> Array which contains the dimensions.</param>
        /// <returns> The equivalent index when the data is expressed as a single array.</returns>
        public static int MultiDimToLinearIdx(int[] indices, int[] Dimensions)
        {
            int idx = 0;

            for (int i = 0; i < Dimensions.Length; i++)
            {
                if (i == Dimensions.Length - 1)
                {
                    idx += indices[Dimensions.Length - 1];
                }
                else
                {
                    idx += indices[i] * DimensionProduct(Dimensions, i + 1);
                }
            }

            return idx;
        }

        /// <summary> Converts a single linear index to multidimensional indices.<para/>
        /// See: https://stackoverflow.com/questions/11316490/convert-a-1d-array-index-to-a-3d-array-index </summary>
        /// <param name="Dimensions"> Array which contains the dimensions.</param>
        /// <param name="idx"> The index when the data is expressed in a single array.</param>
        /// <returns> An array contain the multidimensional indices.</returns>
        public static int[] LinearIdxToMultiDim(int idx, int[] Dimensions)
        {
            int[] indices = new int[Dimensions.Length];

            for (int i = 0; i < Dimensions.Length; i++)
            {
                indices[i] = (idx / DimensionProduct(Dimensions, i + 1)) % Dimensions[i];
            }

            return indices;
        }

        /// <summary> Takes arrays of unique x and y values and computes all of the points on the grid formed by these values.
        /// Similar to MATLAB's meshgrid function.</summary>
        /// <param name="x"> The array of unique x values</param>
        /// <param name="y"> The array of unique y values</param>
        /// <returns>A tuple contained the meshed x and y values represented in row-major order.</returns>
        public static (double[], double[]) MeshGrid1D(double[] X, double[] Y)
        {
            int N = X.Length * Y.Length;

            // Initialise the output
            double[] XX = new double[N];
            double[] YY = new double[N];

            for (int i = 0; i < X.Length; i++)
            {
                for (int j = 0; j < Y.Length; j++)
                {
                    // Find the linear index
                    int idx = j + i * Y.Length;

                    XX[idx] = X[i];
                    YY[idx] = Y[j];
                }
            }

            return (XX, YY);
        }
		
		/// <summary> Returns the mapping from the old indices to the new indices when the data is shifted along a dimension for the FFT.<paramref name="Dim"/>.</summary>
        /// <param name="Indices"> Array of indices arranged in row-major order.</param>
        /// <param name="Dim"> The dimension to apply the shift along.</param>
        /// <param name="Dimensions">Array containing the dimensions of the data.</param>
        /// <returns>The indices with the shift applied.</returns>
        private static int[] GetFFTShiftIndices(int[] Indices, int Dim, int[] Dimensions)
        {
            Indices[Dim] = (Indices[Dim] + (int)Math.Floor(Dimensions[Dim] / 2.0)) % Dimensions[Dim];

            return Indices;
        }

        /// <summary> Returns the mapping from the old indices to the new indices when the data is shifted along a dimension for the IFFT.<paramref name="Dim"/>.</summary>
        /// <param name="Indices"> Array of indices arranged in row-major order.</param>
        /// <param name="Dim"> The dimension to apply the shift along.</param>
        /// <param name="Dimensions">Array containing the dimensions of the data.</param>
        /// <returns>The indices with the shift applied.</returns>
        private static int[] GetIFFTShiftIndices(int[] Indices, int Dim, int[] Dimensions)
        {
            Indices[Dim] = (Indices[Dim] + (int)Math.Ceiling(Dimensions[Dim] / 2.0)) % Dimensions[Dim];

            return Indices;
        }
		
        /// <summary> Shifts the zero-frequency component to the centre of each dimension of the array after applying an IFFT.
        /// See: https://stackoverflow.com/questions/38230794/how-to-write-fftshift-and-ifftshift-in-r 
        /// See: https://au.mathworks.com/help/matlab/ref/fftshift.html 
        /// </summary>
        /// <param name="Input"> Array of multidimensional data arranged in row-major order.</param>
        /// <param name="Dimensions">Array containing the dimensions of the data.</param>
        /// <returns>The data with the zero-frequency shift applied.</returns>
        public static Complex[] IFFTShift(Complex[] Input, int[] Dimensions)
        {
            // Initialise the output
            Complex[] Output = new Complex[Input.Length];
            Complex[] OutputOld = Input.Clone() as Complex[];

            // Loop through the dimensions in reverse
            for (int dim = Dimensions.Length - 1; dim >= 0; dim--)
            {
                // Apply the shift along the current dimension
                for (int i = 0; i < Input.Length; i++)
                {
                    // Convert the 1D index to x, y, z indices
                    int[] indices = LinearIdxToMultiDim(i, Dimensions);

                    // Find the new indicies of the shifted point
                    int[] shifted_indices = GetIFFTShiftIndices(indices, dim, Dimensions);

                    // Find the new 1D index
                    int idx_shifted = MultiDimToLinearIdx(shifted_indices, Dimensions);

                    Output[idx_shifted] = OutputOld[i];
                }

                // Update the old output
                OutputOld = Output.Clone() as Complex[];
            }

            return Output;
        }

        /// <summary> Shifts the zero-frequency component to the centre of each dimension of the array after applying a FFT.
        /// See: https://stackoverflow.com/questions/38230794/how-to-write-fftshift-and-ifftshift-in-r 
        /// See: https://au.mathworks.com/help/matlab/ref/fftshift.html 
        /// </summary>
        /// <param name="Input"> Array of multidimensional data arranged in row-major order.</param>
        /// <param name="Dimensions">Array containing the dimensions of the data.</param>
        /// <returns>The data with the zero-frequency shift applied.</returns>
        public static Complex[] FFTShift(Complex[] Input, int[] Dimensions)
        {
            // Initialise the output
            Complex[] Output = new Complex[Input.Length];
            Complex[] OutputOld = Input.Clone() as Complex[];

            // Loop through the dimensions
            for (int dim = 0; dim < Dimensions.Length; dim++)
            {
                // Apply the shift along the current dimension
                for (int i = 0; i < Input.Length; i++)
                {
                    // Convert the 1D index to x, y, z indices
                    int[] indices = LinearIdxToMultiDim(i, Dimensions);

                    // Find the new indicies of the shifted point
                    int[] shifted_indices = GetFFTShiftIndices(indices, dim, Dimensions);

                    // Find the new 1D index
                    int idx_shifted = MultiDimToLinearIdx(shifted_indices, Dimensions);

                    Output[idx_shifted] = OutputOld[i];
                }

                // Update the old output
                OutputOld = Output.Clone() as Complex[];
            }

            return Output;
        }
		
        /// <summary> Generates the frequency array for a Fourier transform. <para/>
        /// </summary> See: https://dsp.stackexchange.com/questions/7788/setup-frequency-array-properly
        /// <param name="N"> The number of data points in the FFT.</param>
        /// <param name="Fs"> The sampling frequency.</param>
        /// <returns> Array of frequency domain values.</returns>
        public static double[] GetFreqs(int N, double Fs)
        {
            double[] Freqs;

            if (N % 2 == 0)
            {
                Freqs = MathNet.Numerics.Generate.LinearRange(-N / 2, N / 2 - 1);
            }
            else
            {
                Freqs = MathNet.Numerics.Generate.LinearRange(-(N - 1) / 2, (N - 1) / 2);
            }

            for (int i = 0; i < N; i++)
            {
                Freqs[i] *= Fs / N;
            }

            return Freqs;
        }

        /// <summary> Performs cubic interpolation between specified values.<para/>
        /// See: https://stackoverflow.com/questions/20923956/bicubic-interpolation </summary>
        /// <param name="v0"> The first known point.</param>
        /// <param name="v1"> The second known point.</param>
        /// <param name="v2"> The third known point.</param>
        /// <param name="v3"> The fourth known point.</param>
        /// <param name="fracy"> How far between <paramref name="v0"/> and <paramref name="v1"/> to interpolate.
        /// Must be a number between 0 and 1.</param>
        /// <returns> The interpolated value.</returns>
        private static double CubicInterpolate(double v0, double v1, double v2, double v3, double fracy)
        {
            double A = 0.5 * (v3 - v0) + 1.5 * (v1 - v2);
            double B = 0.5 * (v0 + v2) - v1 - A;
            double C = 0.5 * (v2 - v0);
            double D = v1;

            // See Angelo Mascaro's comment for this performance improvement
            return D + fracy * (C + fracy * (B + fracy * A));
        }

        /// <summary> Performs cubic interpolation at a single point.</summary>
        /// <param name="x"> Array of x values.</param>
        /// <param name="F"> Array of function values.</param>
        /// <param name="xval"> x value of the interpolated point.</param>
        /// <returns> The interpolated value.</returns>
        public static double CubicInterpolation(double[] x, double[] F, double xval)
        {
            // Check that x and y have the correct length
            if (x.Length != F.GetLength(0))
            {
                throw new ArgumentException(string.Format("{0} must have the same length as dimension 0 of {1}", nameof(x), nameof(F)), nameof(x));
            }

            // Calculates single point cubic interpolation
            double dx = x[1] - x[0];

            int i = (int)Math.Floor((xval - x[0]) / dx);

            double fracx = (xval - x[i]) / (x[i + 1] - x[i]);

            double p0, p1, p2, p3;

            // Compute the centre data points first
            p1 = F[i]; p2 = F[i + 1];

            // Compute the end points
            if (i == 0)
            {
                p0 = 2 * p1 - p2;
            }
            else
            {
                p0 = F[i - 1];
            }

            if (i == x.Length - 2)
            {
                p3 = 2 * p2 - p1;
            }
            else
            {
                p3 = F[i + 2];
            }

            return CubicInterpolate(p0, p1, p2, p3, fracx);
        }

        /// <summary> Performs bicubic interpolation at a single point. DOES NOT CHECK IF
        /// THE POINT IS WITHIN THE DOMAIN! THIS CHECK WAS REMOVED TO INCREASE PERFORMANCE.</summary>
        /// <param name="x"> Array of x values.</param>
        /// <param name="y"> Array of y values.</param>
        /// <param name="F"> 2D array of function values.</param>
        /// <param name="xval"> x value of the interpolated point.</param>
        /// <param name="yval"> y value of the interpolated point.</param>
        /// <returns> The interpolated value.</returns>
        public static double BicubicInterpolation(double[] x, double[] y, double[,] F, double xval, double yval)
        {
            double dx = x[1] - x[0];
            double dy = y[1] - y[0];

            int i = (int)Math.Floor((xval - x[0]) / dx);
            int j = (int)Math.Floor((yval - y[0]) / dy);

            if (i < 0) { i = 0; } else if (i >= x.Length - 1) { i = x.Length - 2; }
            if (j < 0) { j = 0; } else if (j >= y.Length - 1) { j = y.Length - 2; }

            double fracx = (xval - x[i]) / dx;
            double fracy = (yval - y[j]) / dy;

            double p00, p10, p20, p30;
            double p01, p11, p21, p31;
            double p02, p12, p22, p32;
            double p03, p13, p23, p33;

            /*
            p11 = F[i, j]; p21 = F[i + 1, j];
            p12 = F[i, j + 1]; p22 = F[i + 1, j + 1];
            p01 = F[i - 1, j];
            p02 = F[i - 1, j + 1];
            p31 = F[i + 2, j];
            p32 = F[i + 2, j + 1];
            p10 = F[i, j - 1];
            p20 = F[i + 1, j - 1];
            p13 = F[i, j + 2];
            p23 = F[i + 1, j + 2];
            p00 = F[i - 1, j - 1];
            p03 = F[i - 1, j + 2];
            p30 = F[i + 2, j - 1];
            p33 = F[i + 2, j + 2];
            */

            // Compute the centre cluster first
            p11 = F[i, j]; p21 = F[i + 1, j];
            p12 = F[i, j + 1]; p22 = F[i + 1, j + 1];

            // Compute the centre x edges
            if (i == 0)
            {
                p01 = 2 * p11 - p21;
                p02 = 2 * p12 - p22;
            }
            else
            {
                p01 = F[i - 1, j];
                p02 = F[i - 1, j + 1];
            }

            if (i == x.Length - 2)
            {
                p31 = 2 * p21 - p11;
                p32 = 2 * p22 - p12;
            }
            else
            {
                p31 = F[i + 2, j];
                p32 = F[i + 2, j + 1];
            }

            // Compute the centre y edges
            if (j == 0)
            {
                p10 = 2 * p11 - p12;
                p20 = 2 * p21 - p22;
            }
            else
            {
                p10 = F[i, j - 1];
                p20 = F[i + 1, j - 1];
            }

            if (j == y.Length - 2)
            {
                p13 = 2 * p12 - p11;
                p23 = 2 * p22 - p21;
            }
            else
            {
                p13 = F[i, j + 2];
                p23 = F[i + 1, j + 2];
            }

            // Compute the left corners
            if (i == 0)
            {
                p00 = 2 * p10 - p20;
                p03 = 2 * p13 - p23;
            }
            else
            {
                if (j == 0)
                {
                    p00 = 2 * p10 - p20;
                }
                else
                {
                    p00 = F[i - 1, j - 1];
                }

                if (j == y.Length - 2)
                {
                    p03 = 2 * p13 - p23;
                }
                else
                {
                    p03 = F[i - 1, j + 2];
                }
            }

            // Compute the right corners
            if (i == x.Length - 2)
            {
                p30 = 2 * p20 - p10;
                p33 = 2 * p23 - p13;
            }
            else
            {
                if (j == 0)
                {
                    p30 = 2 * p20 - p10;
                }
                else
                {
                    p30 = F[i + 2, j - 1];
                }

                if (j == y.Length - 2)
                {
                    p33 = 2 * p23 - p13;
                }
                else
                {
                    p33 = F[i + 2, j + 2];
                }
            }

            double x1 = CubicInterpolate(p00, p10, p20, p30, fracx);
            double x2 = CubicInterpolate(p01, p11, p21, p31, fracx);
            double x3 = CubicInterpolate(p02, p12, p22, p32, fracx);
            double x4 = CubicInterpolate(p03, p13, p23, p33, fracx);

            // Calculates single point bicubic interpolation
            return CubicInterpolate(x1, x2, x3, x4, fracy);
        }

        /// <summary> Performs bicubic interpolation for an array of points.</summary>
        /// <param name="x"> Array of x values.</param>
        /// <param name="y"> Array of y values.</param>
        /// <param name="F"> 2D array of function values.</param>
        /// <param name="xvals"> Array of x values for interpolation.</param>
        /// <param name="yvals"> Array of y values for interpolation.</param>
        /// <returns> An array of the interpolated values.</returns>
        public static double[] BicubicInterpolation(double[] x, double[] y, double[,] F, double[] xvals, double[] yvals)
        {
            // Calculates multiple point bilinear interpolation
            double[] zvals = new double[xvals.Length];

            for (int i = 0; i < xvals.Length; i++)
            {
                zvals[i] = BicubicInterpolation(x, y, F, xvals[i], yvals[i]);
            }

            return zvals;
        }

        /// <summary> Reshapes an array in row-major order to a 2D array.</summary>
        /// <param name="input"> The array to be reshaped.</param>
        /// <param name="Dimensions"> Array which contains the dimensions.</param>
        /// <returns> The reshaped array.</returns>
        public static T[,] Reshape1DTo2D<T>(T[] input, int[] Dimensions)
        {
            if (Dimensions.Length != 2)
            {
                throw new ArgumentException(string.Format("The length of {0} must be equal to 2", nameof(Dimensions)), nameof(Dimensions));
            }

            int rowCount = Dimensions[0];
            int colCount = Dimensions[1];

            T[,] output = new T[rowCount, colCount];
            if (rowCount * colCount <= input.Length)
            {
                for (int i = 0; i < rowCount; i++)
                {
                    for (int j = 0; j < colCount; j++)
                    {
                        int idx = MultiDimToLinearIdx(new int[] { i, j }, Dimensions);
                        output[i, j] = input[idx];
                    }
                }
            }
            return output;
        }

        /// <summary> Extension method for <see cref="int"/> which constrains
        /// the value to a particular range.</summary>
        /// <param name="min"> Minimum value.</param>
        /// <param name="max"> Maximum value.</param>
        /// <returns> The constrained value.</returns>
        public static int Clamp(this int value, int min, int max)
        {
            return (value <= min) ? min : (value >= max) ? max : value;
        }

        /// <summary> Extension method for <see cref="double"/> which constrains
        /// the value to a particular range.</summary>
        /// <param name="min"> Minimum value.</param>
        /// <param name="max"> Maximum value.</param>
        /// <returns> The constrained value.</returns>
        public static double Clamp(this double value, int min, int max)
        {
            return (value <= min) ? min : (value >= max) ? max : value;
        }

        /// <summary>
        /// Resize the image to the specified width and height.
        /// Source: https://stackoverflow.com/questions/1922040/how-to-resize-an-image-c-sharp
        /// </summary>
        /// <param name="image">The image to resize.</param>
        /// <param name="width">The width to resize to.</param>
        /// <param name="height">The height to resize to.</param>
        /// <returns>The resized image.</returns>
        public static Bitmap ResizeImage(Image image, int width, int height)
        {
            var destRect = new Rectangle(0, 0, width, height);
            var destImage = new Bitmap(width, height);

            destImage.SetResolution(image.HorizontalResolution, image.VerticalResolution);

            using (var graphics = Graphics.FromImage(destImage))
            {
                graphics.CompositingMode = CompositingMode.SourceCopy;
                graphics.CompositingQuality = CompositingQuality.HighQuality;
                graphics.InterpolationMode = InterpolationMode.HighQualityBicubic;
                graphics.SmoothingMode = SmoothingMode.HighQuality;
                graphics.PixelOffsetMode = PixelOffsetMode.HighQuality;

                using (var wrapMode = new System.Drawing.Imaging.ImageAttributes())
                {
                    wrapMode.SetWrapMode(WrapMode.TileFlipXY);
                    graphics.DrawImage(image, destRect, 0, 0, image.Width, image.Height, GraphicsUnit.Pixel, wrapMode);
                }
            }

            return destImage;
        }
    }
}
