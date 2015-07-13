
#include "spline_curve_fitting.h"
#include "read_write_asc.h"

#include <iostream>

int main(int argc, char *argv[])
{

	char inpf[200],*input; 
	argc--;argv++;					//Skip program name arg

	if (argc<1)
	{
		cout<<"Input file:"<<endl;
		cin>>inpf; 
		input = inpf;
	}
	else input    = argv[0];

	string inFileName( input );
	string outFileName1 = inFileName + "_controls.txt";
	string outFileName2 = inFileName + "_spline.txt";

	std::vector<Vector2d> points;
	CReadWriteAsc::readAsc( inFileName, points );

// 	ClosedCubicBSplineCurve curve(0.002);
// 	SplineCurveFitting scf;
// 	scf.fitAClosedCurve( points, curve, 28, 50, 0.01, 0.001 );
// 	CReadWriteAsc::writeAsc( inFileName, points);
// 	CReadWriteAsc::writeAsc( outFileName1, curve.getControls());
// 	CReadWriteAsc::writeAsc( outFileName2, curve.getSamples() );

	OpenCubicBSplineCurve curve(0.002);
	SplineCurveFitting scf;
	scf.fitAOpenCurve( points, curve, 28, 40, 0.01, 0.001  );
	CReadWriteAsc::writeAsc( inFileName, points);
	CReadWriteAsc::writeAsc( outFileName1, curve.getControls());
	CReadWriteAsc::writeAsc( outFileName2, curve.getSamples() );

}
