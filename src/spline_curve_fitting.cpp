#include "bsplinefitting/spline_curve_fitting.h"
#include <Eigen/SVD>


#include <iostream>


void SplineCurveFitting::initClosedControlPoints(
	const vector<Eigen::Vector2d>& points,
									  vector<Eigen::Vector2d>& controlPs,
									  int controlNum)
{
	const double C_M_PI =3.14159265358979323846;
	// compute the initial 12 control points
	controlPs.clear();

	Eigen::Vector2d min_point = points[0];
	Eigen::Vector2d max_point = points[0];

	for( unsigned int i = 0 ;i!= points.size(); ++i) {
		Eigen::Vector2d v = points[i];
		if( min_point.x() > v.x() )  min_point.x() = v.x();
		if( min_point.y() > v.y() )  min_point.y() = v.y();
		if( max_point.x() < v.x() )  max_point.x() = v.x();
		if( max_point.y() < v.y() )  max_point.y() = v.y();
	}

	Eigen::Vector2d dir = ( max_point-min_point )*0.5;
	Eigen::Vector2d cent = min_point + dir;


	Eigen::Vector2d center = min_point + (max_point - min_point)*0.5;
	double radius = (max_point - min_point).norm()/2.0;

	double delta = 2*C_M_PI/controlNum;

	for(size_t i = 0; i < controlNum; i++)
	{
		Eigen::Vector2d point = center+Eigen::Vector2d(
			radius*std::cos(i*delta),radius*std::sin(i*delta));
		controlPs.push_back(point);
  }
}

double SplineCurveFitting::fitAClosedCurve(
						   const vector<Eigen::Vector2d> &points,
						   ClosedCubicBSplineCurve &curve,
						   int controlNum /* = 28 */,
						   int maxIterNum  /*= 30*/,
						   double alpha /* = 0.002*/,
						   double gama /* = 0.002 */,
						   double eplison /* = 0.0001*/)
{
	controlNum = controlNum/4*4;

	// initialize the cube B-spline
	ClosedCubicBSplineCurve* spline = &curve;
	vector<Eigen::Vector2d> controlPs;
	initClosedControlPoints(points, controlPs, controlNum);
	spline->setNewControl( controlPs);

	// update the control point
	// compute P"(t)
	Eigen::MatrixXd pm = spline->getSIntegralSq();
	Eigen::MatrixXd sm = spline->getFIntegralSq();
	// end test

	// find the foot print, will result in error
	std::vector< std::pair<int,double> > parameters;
	double fsd = spline->findFootPrint( points, parameters);
	int iterNum = 0;
	while( fsd > eplison && iterNum < maxIterNum)
	{
		Eigen::MatrixXd ehm(2*controlNum, 2*controlNum);
		Eigen::VectorXd ehv( 2*controlNum);

		ehm.setZero();
		ehv.setZero();

		for( int i = 0; i!= (int)parameters.size(); ++i)
		{
			spline->getDistance_sd( points[i], parameters[i], ehm, ehv );
		}

		// check if ehm, ehv right
		//solve the function
		Eigen::MatrixXd fm = ehm*0.5 + pm*alpha + sm*gama;
		Eigen::VectorXd ehv2 = ehv*0.5;
		Eigen::JacobiSVD<Eigen::MatrixXd> svd(
			fm, Eigen::ComputeThinU | Eigen::ComputeThinV);
		Eigen::VectorXd resultxy = svd.solve(ehv2);

		// update the curve
		for( int i = 0; i<controlNum; i++)
			controlPs[i] = Eigen::Vector2d( resultxy[i], resultxy[i+controlNum]);
		spline->setNewControl( controlPs );
		++ iterNum;

		fsd = spline->findFootPrint( points, parameters );
	}

	return fsd;
}

void SplineCurveFitting::initOpenControlPoints(
	const vector<Eigen::Vector2d>& points, vector<Eigen::Vector2d>& controlPs,
	int controlNum )
{
	// compute the initial control points
	controlPs.clear();

	double gap = double(points.size())/(controlNum-2);
	controlPs.push_back( points[0] );
	for( int i = 0 ;i!= (controlNum-2); ++i )
	{
		controlPs.push_back( points[std::floor(i*gap)] );
	}
	controlPs.push_back( *points.rbegin());
}

double SplineCurveFitting::fitAOpenCurve( const vector<Eigen::Vector2d>& points,
										 OpenCubicBSplineCurve& curve,
										 int controlNum /*= 28*/, int maxIterNum /*= 30*/,
										 double alpha /*= 0.002*/, double beta /*= 0.005*/,
										 double eplison /*= 0.0001*/ )
{
	// initialize the cube B-spline
	OpenCubicBSplineCurve* spline = &curve;
	vector<Eigen::Vector2d> controlPs;
	initOpenControlPoints(points, controlPs, controlNum);
	spline->setNewControl( controlPs);

	// update the control point
	// compute P"(t)
	Eigen::MatrixXd pm = spline->getSIntegralSq();
	Eigen::MatrixXd sm = spline->getFIntegralSq();
	// end test

	// find the foot print, will result in error
	std::vector< std::pair<int,double> > parameters;
	double fsd = spline->findFootPrint( points, parameters);
	int iterNum = 0;
	while( fsd > eplison && iterNum < maxIterNum)
	{
		Eigen::MatrixXd ehm(2*controlNum, 2*controlNum);
		Eigen::VectorXd ehv( 2*controlNum);

		ehm.setZero();
		ehv.setZero();

		// compute h(D)
		for( int i = 0; i< (int)parameters.size(); i++)
		{
			spline->getDistance_sd( points[i],parameters[i], ehm, ehv );
		}
		// check if ehm, ehv right
		//solve the function
		Eigen::MatrixXd fm = ehm*0.5 + pm*alpha + sm*beta;
		Eigen::VectorXd ehv2 = ehv*0.5;
		Eigen::JacobiSVD<Eigen::MatrixXd> svd(fm,
			Eigen::ComputeThinU | Eigen::ComputeThinV);
		Eigen::VectorXd resultxy = svd.solve(ehv2);

		// update the curve
		for( int i = 0; i<controlNum; i++)
			controlPs[i] = Eigen::Vector2d( resultxy[i], resultxy[i+controlNum]);
		spline->setNewControl( controlPs );
		++ iterNum;

		fsd = spline->findFootPrint( points, parameters );
	}

	return fsd;
}
