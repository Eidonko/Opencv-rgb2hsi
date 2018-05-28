#include <string.h>
#include <iostream>
#include <fstream>
#include <map>
#include <utility>              // std::pair
#include <iomanip>              // std::setprecision
#include <limits>               // std::numeric_limits
#include <vector>
#include <algorithm>    // std::max
#include <queue>

// OpenCV includes
#include <iostream>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui_c.h>

using namespace cv;
using namespace std;


#define GSLOG_USR_ERROR(s) fprintf(stderr, "%s\n", s)

void fromRGB2Mat(Mat& img, unsigned char **r, unsigned char **g, unsigned char **b);
void fromMat2RGB(Mat& img, unsigned char **r, unsigned char **g, unsigned char **b);
bool rowIsAllZeroes(int n, double* buf) ;
bool RGB2HSI_combined(unsigned char r, unsigned char g, unsigned char b, unsigned char& hue, unsigned char& saturation, unsigned char& intensity);
bool RGB2HSI_atan(unsigned char r, unsigned char g, unsigned char b, unsigned char& hue, unsigned char& saturation, unsigned char& intensity);
bool RGB2HSI_acos(unsigned char oR, unsigned char oG, unsigned char oB, unsigned char& hue, unsigned char& saturation, unsigned char& intensity);

void badswitch(int c) {
	fprintf(stderr, "Bad switch: %c\n", c);
	exit(1);
}

bool MiniMax(unsigned char const *r, int rocols, int *rMin, int *rMax)
{
  int min=256;
  int max= -1;
  while (rocols > 0)
  {
	  rocols--;
	  if ( r[rocols] < min ) min = r[rocols];
	  if ( r[rocols] > max ) max = r[rocols];
  }
  *rMin = min, *rMax = max;
  return true;
}

void fromRGB2Mat(Mat& img, unsigned char **r, unsigned char **g, unsigned char **b) {
	for(int y = 0; y < img.rows; y++){
		for(int x = 0; x < img.cols; x++){
			Vec3b& intensity = img.at<Vec3b>(y, x);
			uchar& blu = intensity.val[0];
			uchar& gre = intensity.val[1];
			uchar& red = intensity.val[2];

			blu = b[y][x];
			gre = g[y][x];
			red = r[y][x];
		}
	}
}

void fromMat2RGB(Mat& img, unsigned char **r, unsigned char **g, unsigned char **b) {
	for(int y = 0; y < img.rows; y++){
		for(int x = 0; x < img.cols; x++){
			Vec3b& intensity = img.at<Vec3b>(y, x);
			uchar& blu = intensity.val[0];
			uchar& gre = intensity.val[1];
			uchar& red = intensity.val[2];

			b[y][x] = blu;
			g[y][x] = gre;
			r[y][x] = red;
		}
	}
}

int main(int argc, char *argv[])
{
	bool theSuccess;
	int verbose = 0;
	//unsigned int options;
	bool useArccos = false;

	Mat src, dst;


	//options = 0;

	string oFileName; // output grid
	string iFileName; // input grid
	string tFileName; // clock ticks

	int atan_sing = 0, acos_sing = 0, comb_sing = 0;

	// manage command-line args
	if (argc > 1)
		for (int i = 1; i < argc; i++)
			if (argv[i][0] == '-')
				switch (argv[i][1]) {
				case 'v':	verbose = 1; break;
				case 'o':	if (i >= argc - 1) {
								fprintf(stderr, "Valid switches: -v (verbose), -o outputGrid, -i inputGrid, -c (use arcCosine)\n");
								//badSwitchInRgb2Hsi(argv[i][1]);
								return false;
							}
							oFileName = argv[++i]; break;
				case 'i':	if (i >= argc - 1) {
								fprintf(stderr, "Valid switches: -v (verbose), -o outputGrid, -i inputGrid, -c (use arcCosine)\n");
								//badSwitchInRgb2Hsi(argv[i][1]);
								return false;
							}
							iFileName = argv[++i]; break;
				case 'c':	useArccos = true; break;
				default:	fprintf(stderr, "Unknown switch in component rgb2hsi.\n");
						fprintf(stderr, "Valid switches: -v (verbose), -o outputGrid, -i inputGrid, -c (use arcCosine)\n");
							return false;
				}

	if (oFileName.empty()) {
                fprintf(stderr, "rgb2hsi: option -o outputGrid must be specified.");
                return false;
        }
        if (iFileName.empty()) {
                fprintf(stderr, "rgb2hsi: option -i inputGrid must be specified.");
                return false;
        }

        std::cerr << "Opening input file " << iFileName << '\n';
        src = imread(iFileName, cv::IMREAD_COLOR);
        namedWindow("Input image", CV_WINDOW_AUTOSIZE );
        imshow("Input image", src );

        src.copyTo(dst);

	std::cout << "Image " << iFileName << " consists of " << src.channels() << " channels and "
                << src.cols << "x" << src.rows << " pixels.\n";

	const int rows = src.rows;
        const int cols = src.cols;

	// we need to access the R, G, and B band at the same time => three row buffers are needed
	unsigned char *pr, *pg, *pb;
	unsigned char *ph, *ps, *pi;

	unsigned char **r, **g, **b;
        unsigned char *buffr, *buffg, *buffb;

        int rocol = rows*cols;

        r = new unsigned char *[rows];
        g = new unsigned char *[rows];
        b = new unsigned char *[rows];

        buffr = new unsigned char[ rocol ];
        buffg = new unsigned char[ rocol ];
        buffb = new unsigned char[ rocol ];

        for (int i=0, disp=0; i<rows; i++, disp += cols)
                r[i] = & buffr[disp];
        for (int i=0, disp=0; i<rows; i++, disp += cols)
                g[i] = & buffg[disp];
        for (int i=0, disp=0; i<rows; i++, disp += cols)
                b[i] = & buffb[disp];

        fromMat2RGB(dst, r, g, b);




	atan_sing = acos_sing = comb_sing = 0;

	for (int i = 0; i < rows; i++) {
		unsigned char hue, saturation, intensity;

		pr = r[i], pg = g[i], pb = b[i];

		for (int c = 0; c < cols; c++, pr++, pg++, pb++) {
				theSuccess = RGB2HSI_atan(*pr, *pg, *pb, hue, saturation, intensity);
				if (!theSuccess) {
					atan_sing++;
				}
				theSuccess = RGB2HSI_acos(*pr, *pg, *pb, hue, saturation, intensity);
				if (!theSuccess) {
					acos_sing++;
				}
				theSuccess = RGB2HSI_combined(*pr, *pg, *pb, hue, saturation, intensity);
				if (!theSuccess) {
					comb_sing++;
				}

				if (!theSuccess) {
					*pr = 0, *pg = 0, *pb = 0; // color singularities in Black.  cf. {{1}}
				}
				else {
					*pr = hue;
					*pg = saturation;
					*pb = intensity;
				}
		}

	}

	fprintf(stderr, "Singularities of image %s: atan = %d, acos = %d, combined = %d\n", iFileName.c_str(), atan_sing, acos_sing, comb_sing);

	fromRGB2Mat(dst, r, g, b);
        namedWindow( "HSI image", CV_WINDOW_AUTOSIZE );
        imshow("HSI Normalized image", dst );

        imwrite(oFileName, dst );

        delete r, delete g, delete b;
        delete buffr, delete buffg, delete buffb;

	waitKey(0);
	return 0;
} // EoM.R G B 2 H S I


// Returns true when (hue, saturation, intensity) equivalent of (r,g,b) have been identified.
// Returns false when hue could not be identified. Note that saturation and intensity are identified in this case too.
bool RGB2HSI_acos(unsigned char oR, unsigned char oG, unsigned char oB, unsigned char& hue, unsigned char& saturation, unsigned char& intensity) {
	const double PI = 3.14159265358979323846264338328; // std::atan(1.0)*4.0;
	int r, g, b;
	double h, s, i;
	int m;
	double sqrtop;
	double arcc;


	r = (int)oR, g = (int)oG, b = (int)oB;

	if (r > 255 || g > 255 || b > 255) return false;

	i = (double)(r + g + b) / 3.0;

	m = std::min(r, std::min(g, b));

	s = (i > 0) ? 1.0 - m / i : 0.0;

	saturation = int(round(s*255.0));	// from [0..1] in R  to  [0..255] in N
/*	if (saturation == 255.0) {
		std::cout << "rgb = (" << r << "," << g << "," << b << ")" << std::endl;
		std::cout << "i=" << i << ", m=" << m << "s=" << s << std::endl;
	}
	*/

	intensity = (int) round(i);					// intensity is already in [0..255]

	sqrtop = r*r + g*g + b*b - r*g - r*b - g*b;
	if (sqrtop == 0) {
		// std::cerr << "Hue singularity\n";
		return false;
	}
	arcc = (r - g / 2.0 - b / 2.0) / sqrt(sqrtop);
	if (arcc < -1.0 || arcc > 1.0) {
		// std::cerr << "Hue singularity\n";
		return false;
	}

	arcc = acos(arcc) * 360 / (2 * PI);
	h = (g >= b) ? arcc : 360 - arcc;

	hue = int(round(h / 360.0 * 255.0));		// from [0..360] deg  to  [0..255] in N
	return true;
}

bool RGB2HSI_atan(unsigned char r, unsigned char g, unsigned char b, unsigned char& hue, unsigned char& saturation, unsigned char& intensity) {
	const double PI = 3.14159265358979323846264338328; // std::atan(1.0)*4.0;
	const double sq3 = sqrt(3.0);
	double m;

	if (r > 255 || g > 255 || b > 255) return false;

	m = std::min(r, std::min(g, b));

	double i = (r + g + b) / 3.0;
	intensity = (int)i;
	double s = (intensity > 0) ? 1.0 - m / i : 0.0;
	saturation = int(round(s*255.0));	// from [0..1] in R  to  [0..255] in N


	// Hue:
	// Kender's algorithm: file:///C:/Users/defloriv/Downloads/ADA037772.pdf
	//
	double h;
	if ((r>b) && (g>b)) {
		h = PI / 3.0 + atan(sq3*(g - r) / (r - b + g - b));
	}
	else if (g > r) {
		h = PI + atan(sq3*(b - g) / (b - r + g - r));
	}
	else if (b > g) {
		h = 5.0 * PI / 3.0 + atan(sq3*(r - b) / (r - g + b - g));
	}
	else if (r > b) {
		h = 0.0;
	}
	else
		return false;

	hue = int(round(h * 180.0 / PI));

	//hue = round(hue / 360.0 * 255.0);		// from [0..360] deg  to  [0..255] in N
	return true;
}

bool RGB2HSI_combined(unsigned char r, unsigned char g, unsigned char b, unsigned char& hue, unsigned char& saturation, unsigned char& intensity) {
	const double PI = 3.14159265358979323846264338328; // std::atan(1.0)*4.0;
	const double sq3 = sqrt(3.0);
	double m;


	if (r > 255 || g > 255 || b > 255) return false;

	if (RGB2HSI_atan(r, g, b, hue, saturation, intensity))
		return true;

	if (RGB2HSI_acos(r, g, b, hue, saturation, intensity))
		return true;

	return false;

}

bool rowIsAllZeroes(int n, double* buf) {
	int i;
	for (i = 0; i < n; i++) if (buf[i] != 0) break;
	if (i < n) return false;
	return true;
}

