#include "bsplinefitting/open_cubic_b_spline.h"
#include <fstream>

#include <nabo/nabo.h>

void OpenCubicBSplineCurve::setNewControl(
	const vector<Eigen::Vector2d>& controlPs)
{
	clear();
	controls_ = controlPs;
	int numSeg = nb_segment();

	for( unsigned int i = 0; i<numSeg; i++)
	{
		for( double fj = 0; fj <=1.0f; fj+= interval_)
		{
			Parameter temp(i,fj);
			Eigen::Vector2d p = getPos( temp ) ;
			positions_.push_back(p);
		}
	}
}

Eigen::Vector2d OpenCubicBSplineCurve::getPos( const Parameter& para) const
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
		pm(i,0) = controls_[(ki+i)].x()/6.0;
		pm(i,1) = controls_[(ki+i)].y()/6.0;
	}
	Eigen::MatrixXd rm = tm*cm*pm;

	return Eigen::Vector2d( rm(0,0), rm(0,1));
}

Eigen::Vector2d OpenCubicBSplineCurve::getFirstDiff(const Parameter& para) const
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
		pm(i,0) =  controls_[(ki+i)].x()/ 6.0f;
		pm(i,1) =  controls_[(ki+i)].y() / 6.0f;
	}


	Eigen::MatrixXd rm = tm*cm*pm;

	return Eigen::Vector2d( rm(0,0), rm(0,1));

}

Eigen::Vector2d OpenCubicBSplineCurve::getSecondDiff(
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
		pm(i,0) = controls_[(ki+i)].x()/6.0;
		pm(i,1) =  controls_[(ki+i)].y()/6.0;
	}
	Eigen::MatrixXd rm = tm*cm*pm;

	return Eigen::Vector2d( rm(0,0), rm(0,1));

}

// Refer: http://en.wikipedia.org/wiki/Curvature
double  OpenCubicBSplineCurve::getCurvature(const Parameter& para)  const
{
	Eigen::Vector2d fp = getFirstDiff( para );
	Eigen::Vector2d sp = getSecondDiff( para );

	double kappa = abs( fp.x()*sp.y() - sp.x()*fp.y() );
	kappa = kappa / sqrt( pow( ( fp.x()*fp.x()+fp.y()*fp.y()), 3) );

	return kappa;
}

Eigen::Vector2d OpenCubicBSplineCurve::getTangent( const Parameter &para ) const
{
	Eigen::Vector2d p = getFirstDiff(para);
	return p.normalized();
}

Eigen::Vector2d OpenCubicBSplineCurve::getNormal( const Parameter &para ) const
{
	Eigen::Vector2d v = getTangent( para );
	return Eigen::Vector2d( -v.y(), v.x() );
}


Eigen::Vector2d OpenCubicBSplineCurve::getCurvCenter(
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


double OpenCubicBSplineCurve::findFootPrint(
	const vector<Eigen::Vector2d>& givepoints,
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


OpenCubicBSplineCurve::Parameter OpenCubicBSplineCurve::getPara(
	 int index ) const
{
	int num = (int)( positions_.size()/ nb_segment() );
	int ki = index/num;
	double tf = interval_*( index - ki*num );
	return make_pair( ki, tf);
}

Eigen::VectorXd OpenCubicBSplineCurve::getCoffe( const Parameter& para) const
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
		newv[ (ki+i)] = rv(0,i)/6.0f;
	}
	return newv;
}

bool OpenCubicBSplineCurve::checkSameSide(Eigen::Vector2d p1,
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

bool OpenCubicBSplineCurve::checkInside(Eigen::Vector2d p)
{
	int strip = 0.02/interval_;
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

int OpenCubicBSplineCurve::isLeft( Eigen::Vector2d p0, Eigen::Vector2d p1,
	Eigen::Vector2d p2)
{
	return ( (p1.x() - p0.x()) * (p2.y() - p0.y())
		- (p2.x() - p0.x()) * (p1.y() - p0.y()) );
}

Eigen::MatrixXd OpenCubicBSplineCurve::getSIntegralSq()
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
	unsigned int segNum = nb_segment();
	for( int i = 0; i < segNum; i++)
	{
		for( int j = 0; j < 4; j++)
		{
			for( int n = 0; n < 4; n++)
			{
				int kj = (i+j);
				int kn = (i+n);
				pm(kj,kn) += 2*coffm(j,n);
				pm(controlNum+kj,controlNum+kn) += 2*coffm(j,n);
			}
		}
	}
	return pm;
}

Eigen::MatrixXd OpenCubicBSplineCurve::getFIntegralSq()
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
	int segNum = nb_segment();
	for( int i = 0; i < segNum; i++)
	{
		for( int j = 0; j < 4; j++)
		{
			for( int n = 0; n < 4; n++)
			{
				int kj = ( i+j );
				int kn = (i+n);
				pm(kj,kn) += 2*coffm(j,n);
				pm(controlNum+kj,controlNum+kn) += 2*coffm(j,n);
			}
		}
	}
	return pm;
}

void OpenCubicBSplineCurve::getDistance_sd( const Eigen::Vector2d& point,
	const Parameter& para, Eigen::MatrixXd& ehm, Eigen::VectorXd& ehv )
{
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

	Eigen::MatrixXd leftm = 2*am.transpose()*Nkv*Nkv.transpose()*am;  //8*8
	Eigen::VectorXd rightv= 2*am.transpose()*Nkv*Nkv.transpose()*point;  //8*1

	if( !sign )
	{
		d = -d;
		Eigen::MatrixXd tmpm = 2*am.transpose()*Tkv*Tkv.transpose()*am;
		Eigen::VectorXd tmpv = 2*am.transpose()*Tkv*Tkv.transpose()*point;
		leftm += d/(d-rho)*tmpm;
		rightv += d/(d-rho)*tmpv;
	}

	Eigen::Vector2d oldp = point - neip;
	bool isOuter = false;
	double endpoints_thresh = 1e-2;
	if( ki == 0 && tf <= endpoints_thresh && Tkv.dot( oldp) < 0)
	{
		isOuter = true;
	}
	if( ki == nb_segment()-1 && tf > 1 - endpoints_thresh &&  Tkv.dot(oldp) >0 )
	{
		isOuter = true;
	}
	if( isOuter )
	{
		Eigen::MatrixXd leftm_pd = 2*am.transpose()*am;  //8*8
		Eigen::VectorXd rightv_pd= 2*am.transpose()*point;  //8*1
		double cos_theta = std::abs( oldp.dot(Tkv)/oldp.norm() );
		leftm = cos_theta*leftm_pd + (1-cos_theta)*leftm;
		rightv = cos_theta*rightv_pd + (1-cos_theta)*rightv;
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

int OpenCubicBSplineCurve::local2GlobalIdx( int segId, int localIdx )
{
	int globalIdx = 0;
	if( localIdx < 4)
	{
		globalIdx = segId+localIdx;
	}
	else
	{
		globalIdx = nb_control()+segId+(localIdx-4);
	}
	return globalIdx;
}

void OpenCubicBSplineCurve::getDistance_pd( const Eigen::Vector2d& point,
	const Parameter& para, Eigen::MatrixXd& ehm, Eigen::VectorXd& ehv )
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
	Eigen::MatrixXd rv = tv*cm/6.0;   

	Eigen::MatrixXd am(2,8);
	am.setZero();
	am.block(0,0,1,4) = rv.block(0,0,1,4);
	am.block(1,4,1,4) = rv.block(0,0,1,4);

	Eigen::MatrixXd leftm = 2*am.transpose()*am;
	Eigen::VectorXd rightv= 2*am.transpose()*point;

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
