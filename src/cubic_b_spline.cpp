#include "bsplinefitting/cubic_b_spline.h"
#include <fstream>

#include <nabo/nabo.h>

void ClosedCubicBSplineCurve::setNewControl(
	const vector<Eigen::Vector2d>& controlPs)
{
	clear();
	controls_ = controlPs;

	for( unsigned int i = 0; i<nb_control(); i++)
	{
		for( double fj = 0; fj <=1.0f; fj+= interal_)
		{
			Parameter temp(i,fj);
			Eigen::Vector2d p = getPos( temp ) ;
			positions_.push_back(p);
		}
	}
}

Eigen::Vector2d ClosedCubicBSplineCurve::getPos( const Parameter& para) const
{
	Eigen::MatrixXd cm(4,4);
	cm << -1, 3, -3, 1,
		3, -6, 3, 0,
		-3, 0, 3, 0,
		1, 4, 1, 0;

	double tf = para.second;
	int ki = para.first;

	Eigen::MatrixXd  tm(1,4);
	tm << tf*tf*tf, tf*tf, tf, 1;


	int n = nb_control();
	Eigen::MatrixXd pm(4,2);
	for( int i = 0; i < 4; i++)
	{
		pm(i,0) = controls_[(ki+i)%n].x()/6.0;
		pm(i,1) = controls_[(ki+i)%n].y()/6.0;
	}
	Eigen::MatrixXd rm = tm*cm*pm;

	return Eigen::Vector2d( rm(0,0), rm(0,1));
}

Eigen::Vector2d ClosedCubicBSplineCurve::getFirstDiff(
	const Parameter& para) const
{
	Eigen::MatrixXd cm(4,4);
	cm << -1, 3, -3, 1,
		3, -6, 3, 0,
		-3, 0, 3, 0,
		1, 4, 1, 0;

	double tf = para.second;
	int ki = para.first;

	Eigen::MatrixXd  tm(1,4);
	tm << 3*tf*tf,2*tf, 1, 0;

	int n = nb_control();
	Eigen::MatrixXd pm(4,2);
	for( int i = 0; i < 4; i++)
	{
		pm(i,0) =  controls_[(ki+i)%n].x()/ 6.0f;
		pm(i,1) =  controls_[(ki+i)%n].y() / 6.0f;
	}


	Eigen::MatrixXd rm = tm*cm*pm;

	return Eigen::Vector2d( rm(0,0), rm(0,1));

}


Eigen::Vector2d ClosedCubicBSplineCurve::getSecondDiff(
	const Parameter& para ) const
{
	Eigen::MatrixXd cm(4,4);
	cm << -1, 3, -3, 1,
		3, -6, 3, 0,
		-3, 0, 3, 0,
		1, 4, 1, 0;


	double tf = para.second;
	int ki = para.first;
	Eigen::MatrixXd  tm(1,4);
	tm << 6*tf,2, 0, 0;

	int n = nb_control();
	Eigen::MatrixXd pm(4,2);
	for( int i = 0; i < 4; i++)
	{
		pm(i,0) = controls_[(ki+i)%n].x()/6.0;
		pm(i,1) =  controls_[(ki+i)%n].y()/6.0;
	}
	Eigen::MatrixXd rm = tm*cm*pm;

	return Eigen::Vector2d( rm(0,0), rm(0,1));

}

double  ClosedCubicBSplineCurve::getCurvature(const Parameter& para)  const
{
	Eigen::Vector2d fp = getFirstDiff( para );
	Eigen::Vector2d sp = getSecondDiff( para );

	double kappa = abs( fp.x()*sp.y() - sp.x()*fp.y() );
	kappa = kappa / sqrt( pow( ( fp.x()*fp.x()+fp.y()*fp.y()), 3) );

	return kappa;
}



Eigen::Vector2d ClosedCubicBSplineCurve::getTangent(
	const Parameter &para ) const
{
	Eigen::Vector2d p = getFirstDiff(para);
	return p.normalized();

}

Eigen::Vector2d ClosedCubicBSplineCurve::getNormal(
	const Parameter &para ) const
{
	Eigen::Vector2d v = getTangent( para );
	return Eigen::Vector2d( -v.y(), v.x() );

}


Eigen::Vector2d ClosedCubicBSplineCurve::getCurvCenter(
	const Parameter &para) const
{
	Eigen::Vector2d p = getPos(para);

	Eigen::Vector2d fd = getFirstDiff( para );
	Eigen::Vector2d sd = getSecondDiff( para );

	double p1 = ( fd.x()*fd.x() + fd.y()*fd.y() ) * fd.y();
	double p2 =  sd.y()*fd.x()-sd.x()*fd.y() ;
	double alpha = p.x() - p1/p2;

	double p3 = ( fd.x()*fd.x() + fd.y()*fd.y() ) * fd.x();
	double beta = p.y() + p3/p2;


	return Eigen::Vector2d(alpha ,beta);

}


double ClosedCubicBSplineCurve::findFootPrint(const
	vector<Eigen::Vector2d>& givepoints,
									   vector<Parameter>& footPrints) const
{
	footPrints.clear();
	footPrints.resize( givepoints.size(), Parameter(0,0.0) );

	int iKNei = 1;
	int iDim = 2;
	int iNPts = positions_.size();
	double eps = 0;

	const int kMaxNumNeighbors = 1;
	const int kKdTreeDimension = 2;

	Eigen::MatrixXd data_points(kKdTreeDimension, positions_.size());

	for(size_t idx; idx<positions_.size(); ++idx) {
		data_points.col(idx) = positions_[idx];
	}

	Nabo::NNSearchD* spline_nabo = Nabo::NNSearchD::createKDTreeLinearHeap(
		data_points, kKdTreeDimension);

	const unsigned kSearchOptionFlags =
      Nabo::NNSearchD::ALLOW_SELF_MATCH;
  Eigen::MatrixXi result_indices(kMaxNumNeighbors, 1);
  Eigen::MatrixXd distances(kMaxNumNeighbors, 1);


  Eigen::Vector2d queryPt;
	double squareSum = 0.0;
	for( int i = 0 ;i!= (int)givepoints.size(); ++i) {
		queryPt[0] = givepoints[i].x();
		queryPt[1] = givepoints[i].y();

	  spline_nabo->knn(
	      queryPt, result_indices, distances, kMaxNumNeighbors,
	      eps);
		squareSum += distances(0);
		footPrints[i] =  getPara(result_indices(0)) ;
	}

	return squareSum;

}


ClosedCubicBSplineCurve::Parameter ClosedCubicBSplineCurve::getPara(
	int index ) const
{
	int num = (int)( positions_.size()/controls_.size());
	int ki = index/num;
	double tf = interal_*( index - ki*num );
	return make_pair( ki, tf);
}



Eigen::VectorXd ClosedCubicBSplineCurve::getCoffe( const Parameter& para) const
{
	int ki = para.first;
	double tf = para.second;

	Eigen::Matrix4d cm(4,4);
	cm << -1, 3, -3, 1,
		3, -6, 3, 0,
		-3, 0, 3, 0,
		1, 4, 1, 0;

	Eigen::MatrixXd  tv(1,4);
	tv << tf*tf*tf, tf*tf, tf, 1;

	Eigen::MatrixXd rv = tv*cm;

	Eigen::VectorXd newv(nb_control());
	newv.setZero();
	for( int i = 0; i < 4; i++)
	{
		newv[ (ki+i)%nb_control() ] = rv(0,i)/6.0f;
	}
	return newv;
}

bool ClosedCubicBSplineCurve::checkSameSide(Eigen::Vector2d p1,
	Eigen::Vector2d p2 , Eigen::Vector2d neip)
{
	Eigen::Vector2d v1 = p2 - neip;
	Eigen::Vector2d v2 = p1 - neip;
	bool b = true;

	if( v1.x()*v2.x() + v1.y()*v2.y() < 0)
	{
		b = false;
	}

	return  b;
}

bool ClosedCubicBSplineCurve::checkInside(Eigen::Vector2d p)
{
	int strip = 0.02/interal_;
	int    wn = 0;    // the winding number counter
	// loop through all edges of the polygon
	for (int i=0; i< (int)positions_.size(); i+=strip)
	{
		int j = (i+strip)/(int)positions_.size();
		// edge from V[i] to V[j]
		if (positions_[i].y() <= p.y() ) {
			// start y <= P.y
			if ( positions_[j].y() > p.y() )      // an upward crossing
				if (isLeft( positions_[i], positions_[j], p) > 0)  // P left of edge
					++wn;            // have a valid up intersect
		}
		else {                       // start y > P.y (no test needed)
			if ( positions_[j].y() <= p.y())     // a downward crossing
				if ( isLeft( positions_[i], positions_[j], p) <0 ) // P right of edge
					--wn;            // have a valid down intersect
		}
	}
	if(wn == 0)
		return false;
	else
		return true;

}

int ClosedCubicBSplineCurve::isLeft( Eigen::Vector2d p0,
	Eigen::Vector2d p1, Eigen::Vector2d p2)
{
	return ( (p1.x() - p0.x()) * (p2.y() - p0.y())
		- (p2.x() - p0.x()) * (p1.y() - p0.y()) );
}

Eigen::MatrixXd ClosedCubicBSplineCurve::getSIntegralSq()
{
	// compute P"(t)
	int controlNum  = nb_control();
	Eigen::MatrixXd pm(2*controlNum, 2*controlNum);
	pm.setZero();

	Eigen::Matrix2d tIntergrated;
	tIntergrated << 1/3.0, 1/2.0, 1/2.0, 1.0;
	Eigen::Matrix2d tm;
	tm << 6, 0, 0, 2;


	Eigen::MatrixXd cm(2,4);
	cm << -1, 3, -3, 1,
		3, -6, 3, 0;
	cm = cm/6.0;


	Eigen::Matrix4d coffm = cm.transpose()*tm.transpose()*tIntergrated*tm*cm;
	for( int i = 0; i < controlNum; i++)
	{
		for( int j = 0; j < 4; j++)
		{
			for( int n = 0; n < 4; n++)
			{
				int kj = ( i+j ) % controlNum;
				int kn = (i+n) % controlNum;
				pm(kj,kn) += coffm(j,n);
				pm(controlNum+kj,controlNum+kn) += coffm(j,n);
			}
		}
	}
	return pm;
}

Eigen::MatrixXd ClosedCubicBSplineCurve::getFIntegralSq()
{
	// compute P"(t)
	int controlNum  = nb_control();
	Eigen::MatrixXd pm(2*controlNum, 2*controlNum);
	pm.setZero();


	Eigen::Matrix3d tIntergrated;
	tIntergrated << 1/5.0 , 1/4.0, 1/3.0,
		1/4.0, 1/3.0, 1/2.0,
		1/3.0, 1/2.0, 1/1.0;

	Eigen::MatrixXd cm(3,4);
	cm << -1, 3, -3, 1,
		3, -6, 3, 0,
		-3, 0, 3, 0;
	cm = cm/6.0;

	Eigen::Matrix3d tm;
	tm << 3, 0, 0, 0,2, 0, 0, 0,1;

	Eigen::Matrix4d coffm = cm.transpose()*tm.transpose()*tIntergrated*tm*cm;
	for( int i = 0; i < controlNum; i++)
	{
		for( int j = 0; j < 4; j++)
		{
			for( int n = 0; n < 4; n++)
			{
				int kj = ( i+j ) % controlNum;
				int kn = (i+n) % controlNum;
				pm(kj,kn) += coffm(j,n);
				pm(controlNum+kj,controlNum+kn) += coffm(j,n);
			}
		}
	}
	return pm;

}

void ClosedCubicBSplineCurve::getDistance_sd( const Eigen::Vector2d& point,
	const Parameter& para, Eigen::MatrixXd& ehm, Eigen::VectorXd& ehv )
{
	int controlNum = nb_control();

	int ki = para.first;
	double tf = para.second;

	double kappa = getCurvature( para );
	double rho = 10e+6;
	Eigen::Vector2d neip = getPos( para );
	Eigen::Vector2d Tkv = getTangent( para );
	Eigen::Vector2d Nkv = getNormal( para);
	double d =  ( point - neip ).norm() ;
	Eigen::Vector2d Kv(0.0,0.0);
	bool sign = true;
	if( kappa != 0.0f )
	{
		rho = 1/kappa;
		Kv = getCurvCenter( para );
		double ddd =  ( Kv - neip ).norm() ;
		sign = checkSameSide( Kv, point, neip);
	}

	Eigen::Matrix4d cm(4,4);
	cm << -1, 3, -3, 1,
		3, -6, 3, 0,
		-3, 0, 3, 0,
		1, 4, 1, 0;

	Eigen::MatrixXd  tv(1,4);
	tv << tf*tf*tf, tf*tf, tf, 1;
	Eigen::MatrixXd rv = tv*cm/6.0;   //1*4

	Eigen::MatrixXd am(2,8);
	am.setZero();
	am.block(0,0,1,4) = rv.block(0,0,1,4);
	am.block(1,4,1,4) = rv.block(0,0,1,4);

	Eigen::MatrixXd leftm = 2*am.transpose()*Nkv*Nkv.transpose()*am;
	Eigen::VectorXd rightv= 2*am.transpose()*Nkv*Nkv.transpose()*point;

	if( !sign )
	{
		d = -d;
		Eigen::MatrixXd tmpm = 2*am.transpose()*Tkv*Tkv.transpose()*am;
		Eigen::VectorXd tmpv = 2*am.transpose()*Tkv*Tkv.transpose()*point;
		leftm += d/(d-rho)*tmpm;
		rightv += d/(d-rho)*tmpv;
	}

	for( int iRow = 0 ;iRow!=8; ++iRow )
	{
		int iRowG = local2GlobalIdx(ki,iRow);
		for( int jCol = 0; jCol != 8; ++jCol)
		{
			int jColG = local2GlobalIdx(ki,jCol);
			ehm(iRowG,jColG) += leftm(iRow,jCol);
		}
		ehv(iRowG) += rightv(iRow);
	}
}

int ClosedCubicBSplineCurve::local2GlobalIdx( int segId, int localIdx )
{

	int globalIdx = 0;

	if( localIdx < 4)
	{
		globalIdx = segId+localIdx;
		globalIdx = globalIdx%nb_control();
	}
	else
	{
		globalIdx =segId+(localIdx-4);
		globalIdx = globalIdx%nb_control();
		globalIdx += nb_control();
	}
	return globalIdx;
}
