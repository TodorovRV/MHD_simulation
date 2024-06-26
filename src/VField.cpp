#include <cmath>
#include <Eigen/Eigen>
#include <utility>
#include <utils.h>
#include "VField.h"
#include "BField.h"

using Eigen::Vector3d;


VField::VField(Geometry* geometry): r1_(0.0), r2_(0.0), K1_(1.0), K2_(1.0), b1_(1.0), b2_(1.0), b3_(1.0), is_profile_set(false) {
    geometry_ = geometry;
}

void VField::set_profile(double r1, double r2, double K1, double K2, double b1, double b2, double b3) {
    r1_ = r1;
    r2_ = r2;
    K1_ = K1;
    K2_ = K2;
    b1_ = b1;
    b2_ = b2;
    b3_ = b3;
    is_profile_set = true;
}


ApproximateMHDVfield::ApproximateMHDVfield(Geometry* geometry, double betac_phi) : VField(geometry), betac_phi_(betac_phi) {};


Vector3d ApproximateMHDVfield::vf(const Vector3d &point) const {
    double x = point[0];
    double y = point[1];
    double z = point[2];
    double phi = atan2(y, x);
    double v_phi = c*betac_phi_;
    double r_border = geometry_->radius_at_given_distance(point);
    double r_point = sqrt(x*x + y*y);

    double gamma = density_profile(r_point/r_border, r1_, r2_, K1_, K2_, b1_, b2_, b3_);

    if(z > 0) {
        // Dividing by gamma x,y components accounts for relativistic addition of perpendicular velocities
        return {-sin(phi)*v_phi/gamma, cos(phi)*v_phi/gamma, c*sqrt(1. - 1./(gamma*gamma))};
    } else {
        return {sin(phi)*v_phi/gamma, -cos(phi)*v_phi/gamma, -c*sqrt(1. - 1./(gamma*gamma))};
    }
}


BeskinVfield::BeskinVfield(Geometry* geometry, BeskinMHDBField* bfield, double Psi_tot, double sigma_m, double Omega, double theta) : 
    VField(geometry), bfield_(bfield), Psi_tot_(Psi_tot), sigma_m_(sigma_m), Omega_(Omega), theta_(theta), v_phi_factor_(1.0) {};


Vector3d BeskinVfield::vf(const Vector3d &point) const {
    // Calculating Vfield as drift velosity according to Bfield from Beskin et al 2023
    double x = point[0];
    double y = point[1];
    double z = point[2];
    double phi = atan2(y,x);
    double omega = sqrt(x*x+y*y);

    if (omega < c/Omega_){
        return Vector3d(0, 0, 0);
    }

    auto b = bfield_->bf(point);
    double B_x = b[0];
    double B_y = b[1];
    double B_z = b[2];

    double B_omega = B_x*cos(phi)+B_y*sin(phi);
    double B_phi = B_y*cos(phi)-B_x*sin(phi);
    double B2 = B_omega*B_omega+B_phi*B_phi+B_z*B_z;


    auto ef_ = bfield_->ef(point);
    double E_x = ef_[0];
    double E_y = ef_[1];
    double E_z = ef_[2];

    double E_omega = E_x*cos(phi)+E_y*sin(phi);

    double V_z = (c/B2*E_omega*B_phi);
    double V_omega = (-c/B2*E_z*B_phi);
    double V_phi = v_phi_factor_*(-c/B2*(E_z*B_omega-E_omega*B_z));

    // V_phi = 0.1*c;
    if (V_phi*V_phi+V_omega*V_omega+V_z*V_z > c*c){
        V_phi = sqrt(c*c-V_omega*V_omega-V_z*V_z);
    }

    if (z > 0) {
        return V_omega*Vector3d(cos(phi), sin(phi), 0) + V_phi*Vector3d(-sin(phi), cos(phi), 0) + V_z*Vector3d(0, 0, 1.);
    } else {
        return V_omega*Vector3d(cos(phi), sin(phi), 0) + V_phi*Vector3d(-sin(phi), cos(phi), 0) + V_z*Vector3d(0, 0, -1.);
    }
}

void BeskinVfield::set_v_phi_factor(double v_phi_factor){
    v_phi_factor_ = v_phi_factor;
}


ConstFlatVField::ConstFlatVField(double gamma, Geometry* geometry, double betac_phi) :
    VField(geometry), gamma_(gamma), betac_phi_(betac_phi) {};

Vector3d ConstFlatVField::vf(const Vector3d &point) const {
    double x = point[0];
    double y = point[1];
    double z = point[2];
    double phi = atan2(y, x);
    double v_phi = c*betac_phi_;
    double r_border = geometry_->radius_at_given_distance(point);
    double r_point = sqrt(x*x + y*y);

    double gamma = gamma_;
    if(is_profile_set) {
        gamma *= density_profile(r_point/r_border, r1_, r2_, K1_, K2_, b1_, b2_, b3_);
    }

    if(z > 0) {
        // Relativistic addition of velocities Gamma (along jet) and phi-component (perpendicular jet)
        //return {-sin(phi)*v_phi/gamma, cos(phi)*v_phi/gamma, c*sqrt(1. - 1./(gamma*gamma))};
        return {sin(phi)*v_phi/gamma, -cos(phi)*v_phi/gamma, c*sqrt(1. - 1./(gamma*gamma))};
    } else {
        return {sin(phi)*v_phi/gamma, -cos(phi)*v_phi/gamma, -c*sqrt(1. - 1./(gamma*gamma))};
    }
}


ConstCentralVField::ConstCentralVField(double gamma, Geometry* geometry, double betac_phi, Vector3d origin) :
    VField(geometry), gamma_(gamma), betac_phi_(betac_phi), origin_(std::move(origin)) {}

Vector3d ConstCentralVField::vf(const Vector3d &point) const {
    double x = point[0];
    double y = point[1];
    double z = point[2];
    double phi = atan2(y, x);
    double v_phi = c*betac_phi_;
    double r_border = geometry_->radius_at_given_distance(point);
    double r = (point - origin_).norm();

    double gamma = gamma_;
    if(is_profile_set) {
        gamma *= density_profile(r/r_border, r1_, r2_, K1_, K2_, b1_, b2_, b3_);
    }

    double v_r = c*sqrt(1. - 1./(gamma*gamma));
    if(z > 0) {
        return {v_r*(x-origin_[0])/r - sin(phi)*v_phi, v_r*(y-origin_[1])/r + cos(phi)*v_phi, v_r*(z-origin_[2])/r};
    } else {
        return {v_r*(x-origin_[0])/r + sin(phi)*v_phi, v_r*(y-origin_[1])/r - cos(phi)*v_phi, v_r*(z-origin_[2])/r};
    }
};


//ShearedFlatVField::ShearedFlatVField(double gamma_axis, double gamma_border, double r) :
//    gamma_axis_(gamma_axis), gamma_border_(gamma_border), r_(r) {}
//
//Vector3d ShearedFlatVField::vf(const Vector3d &point) const {
//    double x = point[0];
//    double y = point[1];
//    double z = point[2];
//    double r = sqrt(x*x + y*y);
//    double gamma = gamma_axis_-(gamma_axis_-gamma_border_)*r/r_;
//    if(z > 0) {
//        return {0, 0, c*sqrt(1. - 1./(gamma*gamma))};
//    } else {
//        return {0, 0, -c*sqrt(1. - 1./(gamma*gamma))};
//    }
//}

SheathFlatVField::SheathFlatVField(double gamma_spine_0, double gamma_sheath_0,
                                   Geometry* geometry_in, Geometry* geometry_out,
                                   double gamma_spine_v, double gamma_sheath_v,
                                   double betac_phi_in,  double betac_phi_out) :
    VField(geometry_out),
    gamma_spine_0_(gamma_spine_0),
    gamma_spine_v_(gamma_spine_v),
    gamma_sheath_0_(gamma_sheath_0),
    gamma_sheath_v_(gamma_sheath_v),
    geometry_in_(geometry_in),
    betac_phi_in_(betac_phi_in),
    betac_phi_out_(betac_phi_out)  {}

Vector3d SheathFlatVField::vf(const Vector3d &point) const {
    double x = point[0];
    double y = point[1];
    double z = point[2];
    double phi = atan2(y, x);
    double r_border_out = geometry_->radius_at_given_distance(point);
    double r_border_in = geometry_in_->radius_at_given_distance(point);
    double r_point = sqrt(x*x + y*y);

    double gamma;
    double v_phi;
    if (r_point < r_border_in) {
        gamma = gamma_spine_0_ + gamma_spine_v_*abs(z/pc);
        v_phi = c*betac_phi_in_;
    }
    else if (r_point < r_border_out){
        gamma = gamma_sheath_0_ + gamma_sheath_v_*abs(z/pc);
        v_phi = c*betac_phi_out_;
    } else {
        gamma = 1.0;
        v_phi = 0.0;
    }

    if(z > 0) {
        return {-sin(phi)*v_phi/gamma, cos(phi)*v_phi/gamma, c*sqrt(1. - 1./(gamma*gamma))};
        //return {sin(phi)*v_phi, -cos(phi)*v_phi, c*sqrt(1. - 1./(gamma*gamma))};
    } else {
        return {sin(phi)*v_phi/gamma, -cos(phi)*v_phi/gamma, -c*sqrt(1. - 1./(gamma*gamma))};
    }
}


//ShearedCentralVField::ShearedCentralVField(double gamma_axis, double gamma_border, double theta, Vector3d origin) :
//    gamma_axis_(gamma_axis), gamma_border_(gamma_border), theta_(theta), origin_(std::move(origin)) {}
//
//Vector3d ShearedCentralVField::vf(const Vector3d &point) const {
//    double x = point[0];
//    double y = point[1];
//    double z = point[2];
//    double r = (point - origin_).norm();
//    double theta = acos((z-origin_[2])/r);
//    double gamma = gamma_axis_+(gamma_border_-gamma_axis_)*theta/theta_;
//    double v_r = c*sqrt(1. - 1./(gamma*gamma));
//    return {v_r*(x-origin_[0])/r, v_r*(y-origin_[1])/r, v_r*(z-origin_[2])/r};
//}
//
//SheathCentralVField::SheathCentralVField(double gamma_spine, double gamma_sheath, double theta_sheath, Vector3d origin) :
//    gamma_spine_(gamma_spine), gamma_sheath_(gamma_sheath), theta_sheath_(theta_sheath), origin_(std::move(origin)) {}
//
//Vector3d SheathCentralVField::vf(const Vector3d &point) const {
//    double x = point[0];
//    double y = point[1];
//    double z = point[2];
//    double r = (point - origin_).norm();
//    double theta = acos((z-origin_[2])/r);
//    double gamma;
//    if (theta < theta_sheath_) {
//        gamma = gamma_spine_;
//    } else {
//        gamma = gamma_sheath_;
//    }
//    double v_r = c*sqrt(1. - 1./(gamma*gamma));
//    return {v_r*(x-origin_[0])/r, v_r*(y-origin_[1])/r, v_r*(z-origin_[2])/r};
//}