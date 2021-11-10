using System;
using System.Drawing;
using System.Linq;
using System.Numerics;

namespace DiffractionMaskSimulator
{
    /// <summary>
    /// A <c>static</c> class which provides a range of methods for performing
    /// various physical optics calculations numerically.
    /// </summary>
    public static class PhysicalOptics
    {
        // Constants
        /// <summary> Specifies the options for how all Fourier transforms are calculated using <c>MathNet.Numerics</c>.</summary>
        public const MathNet.Numerics.IntegralTransforms.FourierOptions FOptions = MathNet.Numerics.IntegralTransforms.FourierOptions.Matlab; // Setup the fourier transform to be the same as MATLAB

        public static readonly double[] D65Freqs = new double[] { 380,381,382,383,384,385,386,387,388,389,390,391,392,393,394,395,396,397,398,399,400,401,402,403,404,405,406,407,408,409,410,411,
            412,413,414,415,416,417,418,419,420,421,422,423,424,425,426,427,428,429,430,431,432,433,434,435,436,437,438,439,440,441,442,443,444,445,446,447,448,449,450,451,452,453,454,455,
            456,457,458,459,460,461,462,463,464,465,466,467,468,469,470,471,472,473,474,475,476,477,478,479,480,481,482,483,484,485,486,487,488,489,490,491,492,493,494,495,496,497,498,499,
            500,501,502,503,504,505,506,507,508,509,510,511,512,513,514,515,516,517,518,519,520,521,522,523,524,525,526,527,528,529,530,531,532,533,534,535,536,537,538,539,540,541,542,543,
            544,545,546,547,548,549,550,551,552,553,554,555,556,557,558,559,560,561,562,563,564,565,566,567,568,569,570,571,572,573,574,575,576,577,578,579,580,581,582,583,584,585,586,587,
            588,589,590,591,592,593,594,595,596,597,598,599,600,601,602,603,604,605,606,607,608,609,610,611,612,613,614,615,616,617,618,619,620,621,622,623,624,625,626,627,628,629,630,631,
            632,633,634,635,636,637,638,639,640,641,642,643,644,645,646,647,648,649,650,651,652,653,654,655,656,657,658,659,660,661,662,663,664,665,666,667,668,669,670,671,672,673,674,675,
            676,677,678,679,680,681,682,683,684,685,686,687,688,689,690,691,692,693,694,695,696,697,698,699,700,701,702,703,704,705,706,707,708,709,710,711,712,713,714,715,716,717,718,719,
            720,721,722,723,724,725,726,727,728,729,730,731,732,733,734,735,736,737,738,739,740,741,742,743,744,745,746,747,748,749,750,751,752,753,754,755,756,757,758,759,760,761,762,763,
            764,765,766,767,768,769,770,771,772,773,774,775,776,777,778, 779 };

        public static readonly double[] D65WeightsNormalised = new double[] { 0.424197792941600,0.428163486889352,0.432129180837105,0.436094874784857,0.440063003167199,0.444028697114952,0.447994391062704,
            0.451960085010456,0.455928213392798,0.459893907340551,0.463859601288303,0.487717060274166,0.511574519260029,0.535431978245893,0.559289437231756,0.583144461783029,0.607004355203482,
            0.630861814189345,0.654719273175209,0.678574297726482,0.702431756712345,0.709842175605505,0.717255028933255,0.724665447826415,0.732075866719575,0.739486285612735,0.746899138940485,
            0.754309557833645,0.761719976726805,0.769132830054556,0.776543248947716,0.778193795600003,0.779846776686881,0.781497323339168,0.783150304426046,0.784800851078333,0.786453832165210,
            0.788104378817498,0.789754925469785,0.791407906556663,0.793058453208950,0.787330228617752,0.781602004026555,0.775871345000767,0.770143120409569,0.764412461383781,0.758684236792584,
            0.752956012201386,0.747227787610189,0.741497128584401,0.735768903993203,0.751203219296302,0.766635100164811,0.782069415467910,0.797503730771010,0.812935611639519,0.828369926942618,
            0.843804242245717,0.859233688679636,0.874672872851916,0.890104753720425,0.900409715341563,0.910721980266473,0.921026941887612,0.931331903508750,0.941636865129889,0.951949130054799,
            0.962254091675938,0.972559053297077,0.982871318221986,0.993176279843125,0.993855487093845,0.994541997648336,0.995221204899056,0.995910149888138,0.996589357138858,0.997268564389578,
            0.997955074944069,0.998634282194789,0.999320792749280,1,0.997497401241075,0.994992368047559,0.992489769288634,0.989984736095118,0.987472399597831,0.984967366404316,0.982464767645390,
            0.979959734451875,0.977457135692950,0.974952102499434,0.975852843297880,0.976751149661735,0.977659193763952,0.978559934562398,0.979460675360844,0.980358981724699,0.981259722523145,
            0.982167766625362,0.983066072989218,0.983966813787664,0.977931850438076,0.971896887088489,0.965852186000540,0.959817222650953,0.953782259301366,0.947747295951779,0.941712332602191,
            0.935670065948833,0.929635102599246,0.923600139249658,0.924057812952653,0.924525224394008,0.924982898097002,0.925440571799997,0.925900679937581,0.926365656944346,0.926825765081931,
            0.927283438784925,0.927750850226281,0.928208523929275,0.926891494815872,0.925576900137059,0.924252567719884,0.922937973041071,0.921620943927668,0.920306349248855,0.918989320135452,
            0.917664987718277,0.916350393039465,0.915035798360652,0.912479642040738,0.909925920155414,0.907362460531729,0.904806304211815,0.902252582326492,0.899696426006578,0.897142704121254,
            0.894579244497569,0.892023088177655,0.889469366292332,0.891930579663220,0.894391793034109,0.896853006404997,0.899314219775886,0.901768129843003,0.904229343213892,0.906690556584780,
            0.909151769955669,0.911612983326557,0.914076631132036,0.911291637960625,0.908499341485443,0.905714348314032,0.902922051838850,0.900137058667439,0.897354499930619,0.894562203455436,
            0.891777210284025,0.888984913808843,0.886199920637432,0.885895616313633,0.885588877555243,0.885284573231444,0.884977834473054,0.884673530149255,0.884376529129227,0.884069790370837,
            0.883765486047038,0.883458747288648,0.883154442964849,0.879717021323213,0.876286902985347,0.872849481343710,0.869419363005845,0.865981941364208,0.862544519722572,0.859114401384706,
            0.855676979743070,0.852249295839795,0.848809439763568,0.845698232357044,0.842587024950520,0.839475817543996,0.836364610137472,0.833253402730949,0.830142195324425,0.827028553483311,
            0.823917346076787,0.820806138670263,0.817694931263739,0.817232388691564,0.816767411684799,0.816304869112624,0.815839892105859,0.815377349533684,0.814912372526919,0.814449829954744,
            0.813984852947979,0.813522310375804,0.813059767803629,0.807029673323222,0.801002013277406,0.794971918797000,0.788944258751184,0.782916598705368,0.776886504224961,0.770858844179145,
            0.764831184133329,0.758801089652923,0.752773429607107,0.753893269518688,0.755015543864860,0.756135383776441,0.757255223688022,0.758377498034194,0.759499772380366,0.760619612291947,
            0.761739452203529,0.762861726549700,0.763981566461282,0.763635876749445,0.763290187037610,0.762946931760364,0.762601242048528,0.762253117902102,0.761907428190266,0.761564172913020,
            0.761218483201184,0.760872793489348,0.760527103777512,0.758913073644081,0.757301477945240,0.755687447811808,0.754073417678377,0.752461821979536,0.750847791846105,0.749233761712673,
            0.747622166013833,0.746008135880401,0.744396540181560,0.740652379781534,0.736908219381508,0.733166493416072,0.729422333016045,0.725678172616019,0.721936446650583,0.718192286250557,
            0.714448125850530,0.710706399885095,0.706962239485068,0.707310363631495,0.707658487777921,0.708009046358938,0.708357170505364,0.708705294651791,0.709053418798217,0.709401542944643,
            0.709749667091070,0.710097791237496,0.710448349818513,0.707329839108218,0.704213762832513,0.701095252122218,0.697979175846514,0.694860665136219,0.691744588860514,0.688628512584810,
            0.685510001874515,0.682393925598810,0.679275414888515,0.679436087571481,0.679594325819857,0.679754998502823,0.679913236751198,0.680073909434164,0.680232147682540,0.680392820365506,
            0.680551058613882,0.680709296862257,0.680869969545223,0.682620328015717,0.684373120920801,0.686123479391294,0.687876272296378,0.689626630766871,0.691376989237365,0.693129782142448,
            0.694880140612942,0.696630499083435,0.698383291988519,0.694992124604100,0.691603391654271,0.688212224269852,0.684823491320023,0.681434758370195,0.678043590985776,0.674654858035947,
            0.671263690651528,0.667874957701699,0.664483790317280,0.657217003064953,0.649947781378036,0.642678559691119,0.635411772438792,0.628142550751875,0.620875763499548,0.613606541812631,
            0.606337320125714,0.599070532873387,0.591801311186470,0.593403169146950,0.595007461542020,0.596609319502499,0.598211177462978,0.599813035423458,0.601417327818527,0.603019185779007,
            0.604621043739486,0.606222901699966,0.607824759660445,0.610152079128862,0.612476964162688,0.614801849196515,0.617129168664932,0.619454053698758,0.621778938732585,0.624103823766411,
            0.626431143234828,0.628756028268654,0.631083347737071,0.620264720417360,0.609446093097648,0.598627465777936,0.587808838458224,0.576992645573102,0.566174018253391,0.555355390933679,
            0.544536763613967,0.533720570728845,0.522901943409133,0.529930156071602,0.536960803168660,0.543989015831128,0.551019662928187,0.558047875590655,0.565078522687713,0.572106735350181,
            0.579137382447240,0.586165595109708,0.593196242206766,0.597609872119151,0.602025936466126,0.606442000813101,0.610855630725486,0.615271695072461,0.619685324984846,0.624101389331821,
            0.628517453678796,0.632931083591180,0.637347147938156,0.627589934099856,0.617832720261556,0.608077940857846,0.598320727019546,0.588563513181246,0.578806299342946,0.569051519939237,
            0.559294306100936,0.549537092262637,0.539782312858927,0.525202484097056,0.510625089769776,0.496047695442495,0.481470301115214,0.466890472353344,0.452313078026063,0.437735683698782,
            0.423158289371502,0.408580895044221,0.394003500716941,0.411307461785463,0.428611422853985,0.445917818357097,0.463221779425619,0.480528174928732,0.497832135997254,0.515136097065776,
            0.532442492568888,0.549746453637410,0.567050414705932,0.564146134239592,0.561241853773252,0.558335138872321,0.555430858405981,0.552526577939641,0.549619863038710,0.546715582572370,
            0.543808867671439,0.540904587205099 };

        /// <summary>
        /// Apply a thin lens transformation to an electric field. The lens is a circular lens centred at the origin.
        /// <paramref name="x"/>, <paramref name="y"/>, <paramref name="lambda"/>, and <paramref name="FL"/> must all be in the same units.
        /// Source: https://en.wikipedia.org/wiki/Thin_lens
        /// </summary>
        /// <param name="x">The x coordinate of the E field.</param>
        /// <param name="y">The y coordinate of the E field.</param>
        /// <param name="E">The complex electric field.</param>
        /// <param name="lambda">The wavelength of the electromagnetic wave.</param>
        /// <param name="FL">The focal length of the lens.</param>
        /// <returns>The electric field with the lens transformation applied.</returns>
        public static Complex ApplyThinLens(double x, double y, Complex E, double lambda, double FL)
        {
            return E * Complex.Exp(-Complex.ImaginaryOne * Math.PI / (lambda * FL) * (x * x + y * y));
        }

        /// <summary>
        /// Apply a thin lens transformation to an electric field. The lens is a circular lens centred at the origin.
        /// <paramref name="x"/>, <paramref name="y"/>, <paramref name="lambda"/>, and <paramref name="FL"/> must all be in the same units.<para/>
        /// See: https://en.wikipedia.org/wiki/Thin_lens
        /// </summary>
        /// <param name="x">Array of x coordinates of the E field.</param>
        /// <param name="y">Array of y coordinates of the E field.</param>
        /// <param name="E">Array of complex electric field values.</param>
        /// <param name="lambda">The wavelength of the electromagnetic wave.</param>
        /// <param name="FL">The focal length of the lens.</param>
        /// <returns>The electric field with the lens transformation applied.</returns>
        public static Complex[] ApplyThinLensArray(double[] x, double[] y, Complex[] E, double lambda, double FL)
        {
            Complex[] result = new Complex[E.Length];

            for (int i = 0; i < E.Length; i++)
            {
                result[i] = ApplyThinLens(x[i], y[i], E[i], lambda, FL);
            }

            return result;
        }

        /// <summary>
        /// Compute the E field at distance equal to z with the angular spectrum method <para/>
        /// See: https://en.wikipedia.org/wiki/Angular_spectrum_method <para/>
        /// Source: https://github.com/rafael-fuente/Diffraction-Simulations--Angular-Spectrum-Method/blob/main/diffractsim/monochromatic_simulator.py
        /// </summary>
        /// <param name="x">Array of x coordinates of the E field.</param>
        /// <param name="y">Array of y coordinates of the E field.</param>
        /// <param name="E">Array of complex electric field values.</param>
        /// <param name="Dimensions"> Array which contains the dimensions of the E field data.</param>
        /// <param name="lambda">The wavelength of the electromagnetic wave.</param>
        /// <param name="distance">The distance by which the electric field is propagated.</param>
        /// <returns>The propagated electric field.</returns>
        public static (double[], double[], Complex[]) PropagateField(double[] x, double[] y, Complex[] E, int[] Dimensions, double lambda, double distance)
        {
            MathNet.Numerics.IntegralTransforms.Fourier.Forward2D(E, Dimensions[0], Dimensions[1], FOptions);
            Complex[] Spectrum = Util.FFTShift(E, Dimensions);

            double extent_x = x.Max() - x.Min();
            double extent_y = y.Max() - y.Min();

            int Nx = Dimensions[1];
            int Ny = Dimensions[0];

            double Lx = Math.PI * Nx / 2;
            double Ly = Math.PI * Ny / 2;

            double[] kx = MathNet.Numerics.Generate.LinearSpaced(Nx, Math.Floor(-Lx) / (extent_x / 2), Math.Floor(Lx) / (extent_x / 2));
            double[] ky = MathNet.Numerics.Generate.LinearSpaced(Ny, Math.Floor(-Ly) / (extent_y / 2), Math.Floor(Ly) / (extent_y / 2));

            (double[] kXX, double[] kYY) = Util.MeshGrid1D(kx, ky);

            double[] kZZ = new double[kXX.Length];

            double MAG = Math.Pow(2 * Math.PI / lambda, 2);

            for (int i = 0; i < kZZ.Length; i++)
            {
                kZZ[i] = Math.Sqrt(MAG - kXX[i] * kXX[i] - kYY[i] * kYY[i]);
            }

            // Compute the x and y coordinates at the new location
            double[] XX_prop = new double[kXX.Length];
            double[] YY_prop = new double[kYY.Length];

            double ld = lambda * distance;

            for (int i = 0; i < XX_prop.Length; i++)
            {
                XX_prop[i] = kXX[i] * ld;
                YY_prop[i] = kYY[i] * ld;
            }

            // Propagate the angular spectrum by a distance
            Complex[] SpectrumPropagated = new Complex[Spectrum.Length];

            for (int i = 0; i < Spectrum.Length; i++)
            {
                SpectrumPropagated[i] = Spectrum[i] * Complex.Exp(Complex.ImaginaryOne * kZZ[i] * distance);
            }

            Complex[] E_prop = Util.IFFTShift(SpectrumPropagated, Dimensions);
            MathNet.Numerics.IntegralTransforms.Fourier.InverseMultiDim(E_prop, Dimensions, FOptions);

            return (XX_prop, YY_prop, E_prop);
        }

        /// <summary>
        /// Converts a wavelength of light to RGB colour.<para/>
        /// Source: https://academo.org/demos/wavelength-to-colour-relationship/
        /// </summary>
        /// <param name="lambda">The wavelength of light, m.</param>
        /// <param name="MaxIntensity">Maximum magnitude of the colour. Must be between 0 and 1.</param>
        /// <param name="Gamma">Colour gamma value.</param>
        /// <returns>A <see cref="Color"/> object corresponding to the wavelength <paramref name="lambda"/>.</returns>
        public static Color WavelengthToColour(double lambda, double MaxIntensity, double Gamma) 
        {
            // Initialise the output
            double R, G, B;

            // Convert wavelength to nm
            lambda *= 1E9;

            if (lambda >= 380 && lambda < 440)
            {
                R = (440 - lambda) / (440 - 380);
                G = 0.0;
                B = 1.0;
            }
            else if (lambda >= 440 && lambda < 490)
            {
                R = 0.0;
                G = (lambda - 440) / (490 - 440);
                B = 1.0;
            }
            else if (lambda >= 490 && lambda < 510)
            {
                R = 0.0;
                G = 1.0;
                B = (510 - lambda) / (510 - 490);
            }
            else if (lambda >= 510 && lambda < 580)
            {
                R = (lambda - 510) / (580 - 510);
                G = 1.0;
                B = 0.0;
            }
            else if (lambda >= 580 && lambda < 645)
            {
                R = 1.0;
                G = (645 - lambda) / (645 - 580);
                B = 0.0;
            }
            else if (lambda >= 645 && lambda < 781)
            {
                R = 1.0;
                G = 0.0;
                B = 0.0;
            }
            else
            {
                R = 0.0;
                G = 0.0;
                B = 0.0;
            }

            // Let the intensity fall off near the vision limits
            double Rolloff;
            if ((lambda >= 380) && (lambda < 420))
            {
                Rolloff = 0.3 + 0.7 * (lambda - 380) / (420 - 380);
            }
            else if ((lambda >= 420) && (lambda < 701))
            {
                Rolloff = 1.0;
            }
            else if ((lambda >= 701) && (lambda < 781))
            {
                Rolloff = 0.3 + 0.7 * (780 - lambda) / (780 - 700);
            }
            else
            {
                Rolloff = 0.0;
            }

            // Adjust the final RGB values
            R = Math.Round(255.0 * MaxIntensity * Math.Pow(R * Rolloff, Gamma));
            G = Math.Round(255.0 * MaxIntensity * Math.Pow(G * Rolloff, Gamma));
            B = Math.Round(255.0 * MaxIntensity * Math.Pow(B * Rolloff, Gamma));

            return Color.FromArgb((int)R, (int)G, (int)B);
        }

    }
}
