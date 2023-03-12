/*!
* \file example.cpp
*
* \author Ronghai Hu
* Contact: huronghai@ucas.edu.cn
*
* \brief
*		Examples for calculate LAI using path length distribution (GSL version)
*
*		It is convenient to use GNU Scientific Library (GSL) under Linux/Unix.
*		Visual Studio project file and necessary libraries are provided under Windows.
*		Any question, bug report, or request for a no-GSL version, please contact huronghai@ucas.edu.cn.
*		
*	
* \methods
*		Path length distribution method (Hu et al., 2014)
*		
* \Reference:
*		1. Hu, R., Yan, G., Mu, X., & Luo, J. (2014). Indirect measurement of leaf area index
*		on the basis of path length distribution. REMOTE SENSING OF ENVIRONMENT, 155, 239-247.
*		doi: http://dx.doi.org/10.1016/j.rse.2014.08.032
*
*		2. G. Yan, R. Hu, J. Luo, M. Weiss, H. Jiang, X. Mu, D. Xie, W. Zhang et al. (2019).
*		Review of indirect optical measurements of leaf area index: Recent advances, challenges, and perspectives.
*		Agricultural and Forest Meteorology, 265, 390-411
*		doi:http://dx.doi.org/10.1016/j.agrformet.2018.11.033
*
* \Input
*		Gap fraction
*		Zenith angle 
*		Leaf projection function (G) 
*		Path length distribution
*
* \Output
*		LAI
*		Clumping Index
*
*/

#include <cstdio>
#include <string.h>
#include <list>

#include "LAIPath.h"

void usage(bool wait = false)
{
	fprintf(stderr, "Path Length Distribution Method\nContact: Ronghai HU (huronghai@ucas.edu.cn)\n\n \
Reference: \n1. Hu, R. et al. (2014).Indirect Measurement of Leaf Area Index \
on the Basis of Path Length Distribution.REMOTE SENS ENVIRON, 155, 239 - 247.\n\n\
2. Yan, G. et al. (2019). Review of indirect optical measurements of leaf area index: Recent advances, challenges, and perspectives.\
AGR FOREST METEOROL, 265, 390-411.\n\n");
	//fprintf(stderr, "laslib in.las out.las\n");
	fprintf(stderr, "LAIPATH -i in.txt -o out.txt\n");
	fprintf(stderr, "LAIPATH -i in.txt\n");
	fprintf(stderr, "LAIPATH -h\n");
	//fprintf(stderr, "laslib -i in.las -o out.las\n");
	//fprintf(stderr, "laslib -ilas -olas < in.las > out.las\n");
	//fprintf(stderr, "laslib -h\n");
	if (wait)
	{
		fprintf(stderr, "<press ENTER>\n");
		getc(stdin);
	}
	//exit(1);
}

char* output_path(char* new_path, const char* original_path, const char* fname_append, const char* new_ext)
{


	char drive[_MAX_DRIVE];
	char dir[_MAX_DIR];
	char fname[_MAX_FNAME];
	char ext[_MAX_EXT];


	_splitpath_s(original_path, drive, dir, fname, ext);
	strcat_s(fname, _MAX_FNAME, fname_append);

	if (new_ext != 0)
	{
		strcpy_s(ext, _MAX_EXT, new_ext);
	}

	_makepath_s(new_path, _MAX_PATH, drive, dir, fname, ext);
	return new_path;
}


int main(int argc, char *argv[])
{
	char fname_in[_MAX_PATH];
	char fname_out[_MAX_PATH];
	fname_out[0] = '\0';


	errno_t err;
	FILE * fin, *fout;
	char tmpline[1000];
	//char tmp[100];
	//double tmpPathLen;

	if (argc == 1)
	{
		usage();
		fprintf(stderr, "%s is better run in the command line\n", argv[0]);

		fprintf(stderr, "enter input file: "); fgets(fname_in, 256, stdin);
		fname_in[strlen(fname_in) - 1] = '\0';

		fprintf(stderr, "enter output file: "); fgets(fname_out, 256, stdin);
		fname_out[strlen(fname_out) - 1] = '\0';

	}
	else
	{
		usage();
		//lasreadopener.parse(argc, argv);
		//laswriteopener.parse(argc, argv);
	}

	//command line
	for (int i = 1; i < argc; i++)
	{
		if (argv[i][0] == '\0')
		{
			continue;
		}
		else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "-help") == 0)
		{
			usage();
		}
		else if (strcmp(argv[i], "-i") == 0 )
		{
			if ((i + 1) >= argc)
			{
				fprintf(stderr, "ERROR: '%s' needs 1 argument: stop\n", argv[i]);
				return 1;
			}
			strcpy_s(fname_in, argv[i + 1]);
			i += 1;
		}
		else if (strcmp(argv[i], "-o") == 0)
		{
			if ((i + 1) >= argc)
			{
				fprintf(stderr, "ERROR: '%s' needs 1 argument: stop\n", argv[i]);
				return 1;
			}
			strcpy_s(fname_out, argv[i + 1]);
			i += 1;
		}
		else
		{
			fprintf(stderr, "ERROR: cannot understand argument '%s'\n", argv[i]);
			usage();
		}
	}



	err = fopen_s(&fin, fname_in, "r");
	if (err != 0)
	{
		fprintf(stderr, "ERROR: could not open '%s' for reading\n", fname_in);
		getc(stdin);
		return 1;	
	}

	if (fname_out[0] == '\0')
	{
		output_path(fname_out, fname_in, "_out", "txt");
	}


	err = fopen_s(&fout, fname_out, "w");
	if (err != 0)
	{
		fprintf(stderr, "ERROR: could not open '%s' for writing\n", fname_out);
		getc(stdin);
		//return 1;

	}

	fprintf(fout, "Path Length Distribution Method\r\nContact: Ronghai HU (huronghai@ucas.edu.cn)\r\n \
Reference: \r\n1. Hu, R.et al. (2014).Indirect Measurement of Leaf Area Index \
on the Basis of Path Length Distribution.REMOTE SENS ENVIRON, 155, 239 - 247.\r\n \
2. Yan, G. et al. (2019). Review of indirect optical measurements of leaf area index: Recent advances, challenges, and perspectives.\
AGR FOREST METEOROL, 265, 390-411.\r\n");



	/****************Input 1: Gap probability ******************/
	/*	It is better to remove large gaps before future process
	*	A simple threshold (about 10 times leaf width)
	*	works well for most natural vegetation
	*
	*	If not done, the model also works but might underestimates a little,
	*	large gaps can be considered as path length = 0.
	*/

	double gap_fraction_of_large_gaps;		// Gap fraction of large gaps (e.g., gaps size > 10 * leaf width) (=0 if large gaps are not removed)
	double gap_fraction_inside_canopy;		// Gap fraction after removing large gaps 
	double zenith;							// Observing zenith angle (degree)
	double G;								// Leaf projection function G
	int mode;								// mode of path length distribution input

	double total_gap_fraction, LAIe, CI;

	fgets(tmpline, 1000, fin);
	sscanf_s(tmpline, "%lf", &gap_fraction_inside_canopy);

	fgets(tmpline, 1000, fin);
	sscanf_s(tmpline, "%lf", &gap_fraction_of_large_gaps);

	fgets(tmpline, 1000, fin);
	sscanf_s(tmpline, "%lf", &zenith);

	fgets(tmpline, 1000, fin);
	sscanf_s(tmpline, "%lf", &G);

	fgets(tmpline, 1000, fin);
	sscanf_s(tmpline, "%d", &mode);

	
	total_gap_fraction = gap_fraction_of_large_gaps + (1 - gap_fraction_of_large_gaps)* gap_fraction_inside_canopy;
	LAIe = -log(total_gap_fraction) / G * cos(zenith*M_PI / 180);
	double LAIe_inside_canopy = -log(gap_fraction_inside_canopy) / G * cos(zenith*M_PI / 180);
	//double LAIe_remove_large_gaps = LAIe_inside_canopy * (1 - gap_fraction_of_large_gaps) ;

	printf("\n Input:\nZenith angle:\t%.1f\nG:\t\t%.1f\n", zenith, G);
	printf("Gap fraction of large gaps:\t%.4f\n", gap_fraction_of_large_gaps);
	printf("Gap fraction inside canopy:\t%.4f\n", gap_fraction_inside_canopy);
	printf("Effective LAI (LAIe) = %.4f\n", LAIe);

	fprintf(fout,"\r\n Input:\nZenith angle:\t%.1f\nG:\t\t%.1f\r\n", zenith, G);
	fprintf(fout, "Gap fraction of large gaps:\t%.4f\r\n", gap_fraction_of_large_gaps);
	fprintf(fout, "Gap fraction inside canopy:\t%.4f\r\n", gap_fraction_inside_canopy);
	fprintf(fout, "Effective LAI (LAIe) = %.4f\r\n", LAIe);


	/************Input 2: path length distribution **************/
	/*	path length data could be obtained from 
	*	(1) detailed gap data (photos, TRAC,...) (Hu et al., 2014)
	*	(2) LiDAR data
	*	(3) Ellipse section assumption
	*	
	*
	*	Path length data will be normalized to [0,1] by Stat_hist,
	*	thus only relative path length data is needed.
	*	e.g., {1, 2, 3} and {0.1, 0.2, 0.3} make no difference.
	*/

	double LAI_path;
	gsl_histogram * gsl_hist_path;

	if (mode < 0)	//input path lengths
	{
		/************Input format 1: path lengths **************/
		/*	Path length data will be normalized to[0, 1] by Stat_hist,
		*	thus only relative path length data is needed here.
		*	e.g., { 1, 2, 3 } and {0.1, 0.2, 0.3} make no difference.
		*/

		printf("\nInput mode of path length distribution: Path lengths (-1)\n");
		fprintf(fout, "\r\nInput mode of path length distribution: Path lengths (-1)\r\n");

		int num_of_lines = 0;
		while ( fgets(tmpline, 1000, fin) && tmpline[0] != '\n' )
			++num_of_lines;

		double *path_lengths = new double[num_of_lines];				//path lengths
		
		rewind(fin);
		fgets(tmpline, 1000, fin);
		fgets(tmpline, 1000, fin);
		fgets(tmpline, 1000, fin);
		fgets(tmpline, 1000, fin);
		fgets(tmpline, 1000, fin);

		int num_of_path_lengths = 0;
		double tmp_path_length;
		while (fgets(tmpline, 1000, fin) && tmpline[0] != '\n')
		{
			sscanf_s(tmpline, "%lf", &tmp_path_length);			// read path lengths
			if (tmp_path_length < -1e-6)
			{ 
				fprintf(stderr, "Warning: Path length %lf < 0 will be processed as large gaps\n", tmp_path_length);
				fprintf(fout, "Warning: Path length %lf < 0 will be processed as large gaps\r\n", tmp_path_length);
			}
			else
				path_lengths[num_of_path_lengths++] = tmp_path_length;
		}

		fclose(fin);
		fin = 0;

		gsl_hist_path = gsl_histogram_alloc(NUM_BINS);			// path length distribution
		

		Stat_hist(path_lengths, num_of_path_lengths, gsl_hist_path);			// running statistics to get path length distribution
		delete[] path_lengths;

		printf("\nPath length distribution:\nmin  max  probability\n");
		fprintf(fout, "\r\nPath length distribution:\r\nmin  max  probability\r\n");
		gsl_histogram_fprintf(stdout, gsl_hist_path, "%.2f", "%.3f");
		gsl_histogram_fprintf(fout, gsl_hist_path, "%.2f", "%.3f");

		LAI_path = LAI_PATH(gsl_hist_path, gap_fraction_inside_canopy, zenith, G) \
			* (1 - gap_fraction_of_large_gaps) / num_of_lines * num_of_path_lengths;

		//Fix 2020-03-12: fix the too high estimates when too much path lengths close to 0 observed in path length distribution by Ronghai HU
		if (gsl_fcmp(LAI_path, LAI_MAX, 1e-6) == 0)
		{
			gap_fraction_of_large_gaps = gsl_hist_path->bin[0] / NUM_BINS;
			gsl_hist_path->bin[0] = 0.0;
			gsl_histogram_scale(gsl_hist_path, 1 / (1 - gap_fraction_of_large_gaps));
			LAI_path = LAI_PATH(gsl_hist_path, gap_fraction_inside_canopy, zenith, G) \
				* (1 - gap_fraction_of_large_gaps) / num_of_lines * num_of_path_lengths;

		}
		
	}
	else if (mode > 0)	//input path length distributions 
	{
		/************Input format 2: distribution **************/
		/*	"mode" is the number of bins
		*	Frequency data will be normalized by gsl_histogram_scale
		*	only relative data is needed here.
		*	e.g., { 100, 200, 300 } and {0.1, 0.2, 0.3} make no difference.
		*/

		printf("\nInput mode of path length distribution: Distribution (%d bins)\n", mode);
		fprintf(fout, "\r\nInput mode of path length distribution: Distribution (%d bins)\r\n", mode);

		gsl_hist_path = gsl_histogram_alloc(mode);  // path length distribution
		gsl_histogram_set_ranges_uniform(gsl_hist_path, 0, 1);

		int i = 0;
		while (fgets(tmpline, 1000, fin) && tmpline[0] != '\n')
		{
			sscanf_s(tmpline, "%lf", &gsl_hist_path->bin[i++]);
		}

		fclose(fin);
		fin = 0;

		gsl_histogram_scale(gsl_hist_path, static_cast<double>(gsl_hist_path->n) / gsl_histogram_sum(gsl_hist_path));

		printf("\nPath length distribution:\nmin  max  probability\n");
		fprintf(fout, "\r\nPath length distribution:\r\nmin  max  probability\r\n");
		gsl_histogram_fprintf(stdout, gsl_hist_path, "%.2f", "%.3f");
		gsl_histogram_fprintf(fout, gsl_hist_path, "%.2f", "%.3f");

		LAI_path = LAI_PATH(gsl_hist_path, gap_fraction_inside_canopy, zenith, G) \
			* (1 - gap_fraction_of_large_gaps) ;

	} 
	else    //ellipse section assumption
	{
		/************Input format 3: no input, ellipse section assumption **************/
		/*	Attention: only used when path length distribution is unavailable.
		*	i.e., only effective LAI or gap fraction is available (e.g., LAI-2000)
		*
		*/

		fclose(fin);
		fin = 0;

		printf("\nInput mode of path length distribution: no input, ellipse assumption (0)\n");
		fprintf(fout, "\r\nInput mode of path length distribution: no input, ellipse assumption (0)\r\n");

		LAI_path = LAI_PATH_Circle(gap_fraction_inside_canopy, zenith, G)  \
			* (1 - gap_fraction_of_large_gaps);

		double LAI_path_circle_assumption2 = LAIe2LAI_PATH_Circle(LAIe_inside_canopy, zenith, G)	\
			* (1 - gap_fraction_of_large_gaps);

	}

	CI = LAIe / LAI_path;

	printf("\nResult: LAI_PATH = %.2f\t\tClumping Index = %.3f\n\n", LAI_path, CI);
	fprintf(fout, "\r\nResult: LAI_PATH = %.2f\t\tClumping Index = %.3f\r\n\r\n", LAI_path, CI);

	fclose(fout);
	fout = 0;

	//getc(stdin);

	return 0;
}

