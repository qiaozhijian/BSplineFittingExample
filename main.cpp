//
// Created by qzj on 23-1-23.
//
#include "bsplinefitting/open_cubic_b_spline.h"
#include "bsplinefitting/spline_curve_fitting.h"
#include <iostream>

vector<Eigen::Vector2d> gen_sample_pts(OpenCubicBSplineCurve& spline, int num);

int main(int argc, char** argv)
{
    OpenCubicBSplineCurve* spline = new OpenCubicBSplineCurve();
    vector<Eigen::Vector2d> ctrl_pts;
    ctrl_pts.push_back(Eigen::Vector2d(-2, 1));
    ctrl_pts.push_back(Eigen::Vector2d(0, -0.5));
    ctrl_pts.push_back(Eigen::Vector2d(2, 0.5));
    ctrl_pts.push_back(Eigen::Vector2d(1, 3));
    spline->setNewControl(ctrl_pts);
    vector<Eigen::Vector2d> pts = gen_sample_pts(*spline, 10);

    SplineCurveFitting fitting;
    OpenCubicBSplineCurve* fitted_spline = new OpenCubicBSplineCurve();
    fitting.fitAOpenCurve(pts, *fitted_spline, 4);

    vector<Eigen::Vector2d> fitted_ctrl_pts = fitted_spline->getControls();
    std::cout << "fitted control points:" << std::endl;
    for (int i = 0; i < fitted_ctrl_pts.size(); ++i)
    {
        cout << fitted_ctrl_pts[i].transpose() << endl;
    }

    return 0;
}

vector<Eigen::Vector2d> gen_sample_pts(OpenCubicBSplineCurve& spline, int num)
{
    vector<Eigen::Vector2d> pts;
    for (int i = 0; i < num; ++i)
    {
        OpenCubicBSplineCurve::Parameter para(0, i * 1.0 / num);
        pts.push_back(spline.getPos(para));
    }
    return pts;
}