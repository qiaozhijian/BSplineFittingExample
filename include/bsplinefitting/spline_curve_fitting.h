#include "bsplinefitting/cubic_b_spline.h"
#include "bsplinefitting/open_cubic_b_spline.h"

class  SplineCurveFitting
{
public:
	SplineCurveFitting(void){}
	~SplineCurveFitting(void){}

	//////////////////////////////////////////////////////////////////////////
	// Fitting B-Spline Curves to Point Clouds
	// refer : Fitting B-Spline Curve to Point Clouds by Curvature-Based Squared Distance Minimization
	// controlNum: the number of control points
	// alpha:      the coefficient of curvature constraint
	// beta:       the coefficient of curve length constraint
	// maxIterNum: the maximum of iteration
	// eplison:    the threshold of ending iteration
	//////////////////////////////////////////////////////////////////////////
	double fitAClosedCurve( const vector<Eigen::Vector2d>& points,
		ClosedCubicBSplineCurve& curve,
		int controlNum  = 28,
		int maxIterNum = 30,
		double alpha  = 0.00,
		double beta = 0.005,
		double eplison = 0.0001);


	double fitAOpenCurve( const vector<Eigen::Vector2d>& points,
		OpenCubicBSplineCurve& curve,
		int controlNum  = 28,
		int maxIterNum = 30,
		double alpha  = 0.002,
		double beta = 0.005,
		double eplison = 0.0001);

private:
	// initial control points
	void initClosedControlPoints( const vector<Eigen::Vector2d>& points,
		vector<Eigen::Vector2d>& controlPs,
		int controlNum);


	void initOpenControlPoints(const vector<Eigen::Vector2d>& points,
		vector<Eigen::Vector2d>& controlPs,
		int controlNum);
};
