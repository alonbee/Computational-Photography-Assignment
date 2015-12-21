#include "a2.h"
#include <iostream>
#include <math.h>

using namespace std;


// This is a way for you to test your functions. 
// We will only grade the contents of a2.cpp
int main()
{	
	FloatImage wheat("Input/wheat.png");
	FloatImage equal_input("Input/equal_input.png");
	FloatImage over_exposed("Input/over_exposed.png");
	FloatImage under_exposed("Input/under_exposed.png");

	FloatImage Nelson_Chromatic_Progression("Input/Nelson_Chromatic_Progression.png");	// Test for 
	FloatImage Nelson_matched_wheat("Input/Nelson_Chromatic_Progression_histogram_after.png");


	FloatImage linear_wheat = changeGamma(wheat, 2.2f, 1.0f);
	linear_wheat.write("Output/linear_wheat.png");

	FloatImage brighter_wheat = exposure(wheat, 1.25f);
	brighter_wheat.write("Output/brighter_wheat.png");

	FloatImage grey_wheat = color2gray(wheat);
	grey_wheat.write("Output/grey_wheat.png");

	vector<FloatImage> lc_wheat = lumiChromi(wheat);
	if (lc_wheat.size() == 2)
	{
		lc_wheat[0].write("Output/wheat_lumi.png");
		lc_wheat[1].write("Output/wheat_chromi.png");
	}

	FloatImage brighter_contrastier_wheat = brightnessContrastLumi(wheat, 1.5, 2.0f);
	brighter_contrastier_wheat.write("Output/brighter_contrastier_wheat.png");

	FloatImage yuv_wheat = rgb2yuv(wheat);
	FloatImage rgb_wheat = yuv2rgb(yuv_wheat);

	FloatImage wheat_difference = (rgb_wheat - wheat) * 100.0f;
	for (unsigned int i=0; i < wheat_difference.size(); ++i)
	wheat_difference(i) = fabs(wheat_difference(i));
	wheat_difference.write("Output/wheat_difference.png");

	vector<FloatImage> spanish_wheat = spanish(wheat);
	if (spanish_wheat.size() == 2)
	{
		spanish_wheat[0].write("Output/spanish_wheat_lumi.png");
		spanish_wheat[1].write("Output/spanish_wheat_chromi.png");
	}

	vector<FloatImage> spanish_Nelson_Chromatic_Progression = spanish(Nelson_Chromatic_Progression);	// My castle illusion 
	if (spanish_Nelson_Chromatic_Progression.size() == 2)
	{
		spanish_Nelson_Chromatic_Progression[0].write("Output/spanish_Nelson_Chromatic_Progression_lumi.png");
		spanish_Nelson_Chromatic_Progression[1].write("Output/spanish_Nelson_Chromatic_Progression_chromi.png");
	}

	FloatImage grayworld_flower = grayworld(FloatImage("Input/flower.png"));	
	grayworld_flower.write("Output/grayworld_flower.png");

	FloatImage test_white_balance = grayworld(FloatImage("Input/test_white_balance.png"));	
	test_white_balance.write("Output/grayworld_test_white_balance.png");

	FloatImage over_exposed_histogram = visualizeRGBHistogram(Histogram(over_exposed, 0), Histogram(over_exposed, 1), Histogram(over_exposed, 2));
	over_exposed_histogram.write("Output/overexposed_histogram.png");

	FloatImage under_exposed_histogram = visualizeRGBHistogram(Histogram(under_exposed, 0), Histogram(under_exposed, 1), Histogram(under_exposed, 2));
	under_exposed_histogram.write("Output/underexposed_histogram.png");


	FloatImage wheat_histogram = visualizeRGBHistogram(Histogram(wheat, 0), Histogram(wheat, 1), Histogram(wheat, 2));
	wheat_histogram.write("Output/wheat_histogram.png");

	FloatImage wheat_histogram_equal = equalizeRGBHistograms(wheat);
	wheat_histogram_equal.write("Output/equalized_wheat.png");

	FloatImage equal_input_histogram = visualizeRGBHistogram(Histogram(equal_input, 0), Histogram(equal_input, 1), Histogram(equal_input, 2));
	equal_input_histogram.write("Output/equalized_wheat_histogram.png"); // Output the histogram after equlization. 


	 FloatImage Nelson_Chromatic_Progression_histogram = visualizeRGBHistogram(Histogram(Nelson_Chromatic_Progression, 0), Histogram(Nelson_Chromatic_Progression, 1), Histogram(Nelson_Chromatic_Progression, 2));
	 Nelson_Chromatic_Progression_histogram.write("Output/Nelson_histogram.png");

	 FloatImage Nelson_Chromatic_Progression_histogram_after = matchRGBHistograms(Nelson_Chromatic_Progression,wheat);	// Match my picture(Nelson_Chromatic_Progression.png) with wheat.png.
	 Nelson_Chromatic_Progression_histogram_after.write("Output/Nelson_matched_wheat.png");

	 FloatImage Nelson_matched_wheat_histogram = visualizeRGBHistogram(Histogram(Nelson_matched_wheat, 0), Histogram(Nelson_matched_wheat, 1), Histogram(Nelson_matched_wheat, 2));
	 Nelson_matched_wheat_histogram.write("Output/Nelson_matched_wheat_histogram.png");


}
