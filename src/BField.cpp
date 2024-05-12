#include <cmath>
#include <Eigen/Eigen>
#include <BField.h>
#include <utils.h>
#include <boost/math/special_functions/bessel.hpp>

using Eigen::Vector3d;


VectorBField::VectorBField(bool in_plasma_frame, double tangled_fraction, Geometry* geometry_out, Geometry* geometry_in) :
    in_plasma_frame_(in_plasma_frame),
    tangled_fraction_(tangled_fraction) {
    geometry_in_ = geometry_in;
    geometry_out_ = geometry_out;
}

double VectorBField::get_tangled_fraction(const Vector3d &point) const {
    return tangled_fraction_;
}

Vector3d VectorBField::bf(const Vector3d &point) const {
	double r_border;
    double x = point[0];
    double y = point[1];
    double r_point = sqrt(x*x + y*y);
//    std::cout << "r = " << r_point/pc << "\n";
	
	
	
	double factor = 1.0;
	
	if(geometry_out_) {
		// Find radius of outer surface at given point
		r_border = geometry_out_->radius_at_given_distance(point);
		if (r_point > r_border) {
			factor = exp(-pow(r_point - r_border, 2.0)/l_eps_B/l_eps_B);
		}
	}
//	if(geometry_in_) {
//		// Find radius of inner surface at given point
//		r_border = geometry_in_->radius_at_given_distance(point);
//		if (r_point < r_border) {
//			factor = exp(-pow(r_point - r_border, 2.0)/l_eps_B/l_eps_B);
//		}
//	}
	return _bf(point)*factor;
	

// FIXME: Original
//    if(geometry_out_) {
//        // Find radius of outer surface at given point
//        double r_border_out = geometry_out_->radius_at_given_distance(point);
////        std::cout << "r_border = " << r_border_out/pc << "\n";
//        if (r_point >= r_border_out) {
//            return {0.0, 0.0, 0.0};
//        }
//    }
//    if(geometry_in_) {
//        // Find radius of inner surface at given point
//        std::cout << "NEVER =========================================== " << "\n";
//        double r_border_in = geometry_in_->radius_at_given_distance(point);
//        if (r_point < r_border_in) {
//            return {0.0, 0.0, 0.0};
//        }
//    }
//    return _bf(point);
}
/*
Vector3d VectorBField::bf_plasma_frame(const Vector3d &point, Vector3d &v) const {
    Vector3d b = bf(point);
    if (in_plasma_frame_) {
        return b;
    } else {
        return get_B_prime(b, v);
    }
}*/

double VectorBField::bf_tangled_plasma_frame(const Vector3d &point, Vector3d &v) const {
    Vector3d b = bf(point);
    if (tangled_fraction_ > 0.0) {
        if (in_plasma_frame_) {
            return tangled_fraction_*b.norm();
        } else {
            return tangled_fraction_*get_B_prime(b, v).norm();
        }
    } else {
        return 0.0;
    }
}

Vector3d VectorBField::bhat_lab_frame(const Vector3d &point, Vector3d &v) const {
    Vector3d b = bf(point);
    if (in_plasma_frame_) {
        auto b_hat_prime = b.normalized();
        Vector3d beta = v/c;
        return get_B_hat(b_hat_prime, beta);
    } else {
        return b.normalized();
    }
}


ConstCylinderBFieldZ::ConstCylinderBFieldZ(double b_0, double n_b, bool in_plasma_frame, double tangled_fraction, Geometry* geometry_out, Geometry* geometry_in) :
    VectorBField(in_plasma_frame, tangled_fraction, geometry_out, geometry_in), b_0_(b_0), n_b_(n_b) {};

Vector3d ConstCylinderBFieldZ::_bf(const Vector3d &point) const {
    double r = abs(point[2]);
    double z = point[2];
    if(z > 0) {
        return Vector3d(0.0, 0.0, b_0_*pow(r/pc, -n_b_));
    } else {
        return Vector3d(0.0, 0.0, -b_0_*pow(r/pc, -n_b_));
    }
}


//RadialConicalBField::RadialConicalBField(double b_0, double n_b, bool in_plasma_frame, double tangled_fraction) :
//    VectorBField(in_plasma_frame, tangled_fraction), b_0_(b_0), n_b_(n_b) {};
//
//Vector3d RadialConicalBField::bf(const Vector3d &point) const {
//    double r = point.norm();
//    double z = point[2];
//    return Vector3d(b_0_*pow(r/pc, -n_b_)*point[0]/r,
//                    b_0_*pow(r/pc, -n_b_)*point[1]/r,
//                    b_0_*pow(r/pc, -n_b_)*point[2]/r);
//}
//


ToroidalBField::ToroidalBField(double b_0, double n_b, bool in_plasma_frame, double tangled_fraction, Geometry* geometry_out, Geometry* geometry_in) :
    VectorBField(in_plasma_frame, tangled_fraction, geometry_out, geometry_in), b_0_(b_0), n_b_(n_b){};

Vector3d ToroidalBField::_bf(const Vector3d &point) const {
    double x = point[0];
    double y = point[1];
    double z = point[2];
    double fi = atan2(y, x);
    double b = b_0_*pow(z/pc, -n_b_);
    return {-sin(fi)*b, cos(fi)*b, 0};
}

Vector3d ToroidalBField::bf_plasma_frame(const Vector3d &point, Vector3d &v) const {
    Vector3d b = bf(point);
    if (in_plasma_frame_) {
        return b;
    } else {
        return get_B_prime(b, v);
    }
}

// TODO: Check helicity (rotation direction) for z < 0
HelicalCylinderBField::HelicalCylinderBField(double b_0, double pitch_angle, bool in_plasma_frame, double tangled_fraction, Geometry* geometry_out, Geometry* geometry_in) :
    VectorBField(in_plasma_frame, tangled_fraction, geometry_out, geometry_in), b_0_(b_0), pitch_angle_(pitch_angle) {};

Vector3d HelicalCylinderBField::_bf(const Vector3d &point) const {
    double z = point[2];
    double phi = atan2(point[1], point[0]);
    double b_z = b_0_*cos(pitch_angle_);
    // FIXME: Here was + sign
    double b_phi = b_0_*sin(pitch_angle_);
    if(z > 0) {
        return b_phi*Vector3d(-sin(phi), cos(phi), 0) + b_z*Vector3d(0, 0, -1);
    } else {
        return b_phi*Vector3d(-sin(phi), cos(phi), 0) + b_z*Vector3d(0, 0, -1);
    }
}


// TODO: Check helicity (rotation direction) for z < 0
HelicalConicalBField::HelicalConicalBField(double b_0, double n_b, double pitch_angle, bool in_plasma_frame, double tangled_fraction, Geometry* geometry_out, Geometry* geometry_in) :
    VectorBField(in_plasma_frame, tangled_fraction, geometry_out, geometry_in), b_0_(b_0), n_b_(n_b), pitch_angle_(pitch_angle) {};

Vector3d HelicalConicalBField::_bf(const Vector3d &point) const {
    double z = point[2];
    double b = b_0_*pow(z/pc, -n_b_);
    double phi = atan2(point[1], point[0]);
    double b_z = b*cos(pitch_angle_);
    // double b_phi = b*sin(pitch_angle_);
    double b_phi = b*sin(pitch_angle_)*sqrt(point[0]*point[0]+ point[1]*point[1])/geometry_out_->radius_at_given_distance(point);
    if(z > 0) {
        return b_phi*Vector3d(-sin(phi), cos(phi), 0) + b_z*Vector3d(0, 0, 1);
    } else {
        return b_phi*Vector3d(-sin(phi), cos(phi), 0) + b_z*Vector3d(0, 0, -1);
    }
}

Vector3d HelicalConicalBField::bf_plasma_frame(const Vector3d &point, Vector3d &v) const {
    Vector3d b = bf(point);
    if (in_plasma_frame_) {
        return b;
    } else {
        return get_B_prime(b, v);
    }
}


//SpiralConicalBField::SpiralConicalBField(double b_0, double pitch_angle, bool in_plasma_frame, double tangled_fraction) :
//    VectorBField(in_plasma_frame, tangled_fraction), b_0_(b_0), pitch_angle_(pitch_angle) {};
//
//Vector3d SpiralConicalBField::bf(const Vector3d &point) const {
//    double z = point[2];
//    double x = point[0];
//    double y = point[1];
//    double b_z = b_0_/(z*z/(pc*pc));
//    return Vector3d(b_z*(x/z + y*tan(pitch_angle_)/pc),
//                    b_z*(y/z - x*tan(pitch_angle_)/pc),
//                    b_z);
//}
//
//ForceFreeCylindricalBField::ForceFreeCylindricalBField(double b_0, double mu, bool in_plasma_frame, double tangled_fraction) :
//    VectorBField(in_plasma_frame, tangled_fraction), b_0_(b_0), mu_(mu) {};
//
//Vector3d ForceFreeCylindricalBField::bf(const Vector3d &point) const {
//    double x = point[0];
//    double y = point[1];
//    double atan_term = atan(y/x);
//    double bessel_0 = boost::math::cyl_bessel_i(0, mu_);
//    double bessel_1 = boost::math::cyl_bessel_i(1, mu_);
//    return Vector3d(-b_0_*bessel_1*sin(atan_term),
//                    b_0_*bessel_1*sin(atan_term),
//                    b_0_*bessel_0);
//}


ReversedPinchCylindricalBField::ReversedPinchCylindricalBField(double b_0, double tangled_fraction,
                                                               Geometry* geometry_out, Geometry* geometry_in) :
    VectorBField(true, tangled_fraction, geometry_out, geometry_in), b_0_(b_0) {};

Vector3d ReversedPinchCylindricalBField::_bf(const Vector3d &point) const {
    double x = point[0];
    double y = point[1];
    double phi = atan2(y, x);
    double r_border = geometry_out_->radius_at_given_distance(point);
    double ro_normed = sqrt(point[0]*point[0]+ point[1]*point[1])/r_border;

    double bessel_0 = boost::math::cyl_bessel_i(0, 2.405*ro_normed);
    double bessel_1 = boost::math::cyl_bessel_i(1, 2.405*ro_normed);
    auto B_cyl =  Vector3d(0, b_0_*bessel_1,b_0_*bessel_0);

    Eigen::Matrix3d cyl_to_cart;
    cyl_to_cart << cos(phi), -sin(phi), 0,
                   sin(phi), cos(phi), 0,
                   0, 0, 1;
    return cyl_to_cart*B_cyl;
}


ReversedPinchConicalBField::ReversedPinchConicalBField(double b_0, double n_b, Geometry* geometry, double tangled_fraction) :
    VectorBField(true, tangled_fraction), b_0_{b_0}, n_b_(n_b) {
    geometry_ = geometry;
}

Vector3d ReversedPinchConicalBField::_bf(const Vector3d &point) const {
    double x = point[0];
    double y = point[1];
    double z = point[2];
    double b = b_0_*pow(z/pc, -n_b_);
    double phi = atan2(y, x);
    // Find radius at given point
    double r_border = geometry_->radius_at_given_distance(point);
    double ro_normed = sqrt(x*x + y*y)/r_border;

    double bessel_0 = boost::math::cyl_bessel_i(0, 2.405*ro_normed);
    double bessel_1 = boost::math::cyl_bessel_i(1, 2.405*ro_normed);
    auto B_cyl =  Vector3d(0, b*bessel_1,b*bessel_0);

    Eigen::Matrix3d cyl_to_cart;
    cyl_to_cart << cos(phi), -sin(phi), 0,
        sin(phi), cos(phi), 0,
        0, 0, 1;

    return cyl_to_cart*B_cyl;
}


BeskinMHDBField::BeskinMHDBField(double Psi_tot, double sigma_m, double Omega, double theta) : 
    VectorBField(false, 0), Psi_tot_(Psi_tot), sigma_m_(sigma_m), Omega_(Omega), theta_(theta) {};

Vector3d BeskinMHDBField::_bf(const Vector3d &point) const {
    double x = point[0];
    double y = point[1];
    double z = abs(point[2]);
    double phi = atan2(y,x);
    double gamma = 1.;
    double omega = sqrt(x*x+y*y);
    double power = 1.-0.025*log(z*abs(Omega_)*theta_/c/sqrt(gamma*sigma_m_));
    double br1 = 1.+(omega*omega*Omega_*Omega_/c/c/gamma/gamma);
    double br2 = 1.+(z*z*Omega_*Omega_*theta_*theta_/c/c/gamma/sigma_m_);
    double denom = 4.*sigma_m_*(-1.+pow(br2, power))*(-1.+pow(br2, power))*power*power;
    double grad_psi_z = Psi_tot_*gamma*gamma*(-0.05/z*(-1.+pow(br2, power))*(-1.+pow(br1, power))*2.*power+
            0.1/z*(-1.+pow(br2, power))*(-1.+pow(br1, power))*power-
            4.*pow(br2, power)*(-1.+pow(br1, power))*power*power*
            (2.*z*Omega_*Omega_*theta_*theta_*power/(c*c*gamma*sigma_m_+z*z*Omega_*Omega_*theta_*theta_)-0.025*log(br2)/z)-
            0.1/z*pow(br1, power)*(-1.+pow(br2, power))*power*power*log(br1))/denom;
    double grad_psi_omega = Psi_tot_*8.*omega*Omega_*Omega_/c/c*pow(br1, power-1.)*(-1.+pow(br2, power))*power*power*power/denom;
    double psi = gamma*gamma*Psi_tot_*(-1.+pow(br1, power))/sigma_m_/(-1.+pow(br2, power));
    double Gamma_max = gamma + 2.*sigma_m_*psi/Psi_tot_*(1.-psi/Psi_tot_);
    double epsilon = 0.5/Gamma_max/Gamma_max;
    double Omega_F = Omega_*sqrt(1.-psi/Psi_tot_);

    double B_z;
    if (omega == 0){
        B_z = Psi_tot_*8.*Omega_*Omega_/c/c*pow(br1, power-1.)*(-1.+pow(br2, power))*power*power*power/denom/2./pi;
        if (z > 0) {
            return B_z*Vector3d(0, 0, 1.);
        } else {
            return B_z*Vector3d(0, 0, -1.);
        }
    }

    double B_omega = grad_psi_z/2./pi/omega;
    B_z = grad_psi_omega/2./pi/omega;
    double B_phi = -(1.+epsilon)*Omega_F/2./pi/c*sqrt(grad_psi_z*grad_psi_z+grad_psi_omega*grad_psi_omega);
    if (Psi_tot_ < 0){
        B_phi = -B_phi;
    }
    if (point[2] > 0) {
        return B_omega*Vector3d(cos(phi), sin(phi), 0) + B_phi*Vector3d(-sin(phi), cos(phi), 0) + B_z*Vector3d(0, 0, 1.);
    } else {
        return B_omega*Vector3d(cos(phi), sin(phi), 0) + B_phi*Vector3d(-sin(phi), cos(phi), 0) + B_z*Vector3d(0, 0, -1.);
    }
}

Vector3d BeskinMHDBField::ef(const Vector3d &point) const {
    double x = point[0];
    double y = point[1];
    double z = abs(point[2]);
    double phi = atan2(y,x);
    double gamma = 1.;
    double omega = sqrt(x*x+y*y);
    double power = 1.-0.025*log(z*abs(Omega_)*theta_/c/sqrt(gamma*sigma_m_));
    double br1 = 1.+(omega*omega*Omega_*Omega_/c/c/gamma/gamma);
    double br2 = 1.+(z*z*Omega_*Omega_*theta_*theta_/c/c/gamma/sigma_m_);
    double denom = 4.*sigma_m_*(-1.+pow(br2, power))*(-1.+pow(br2, power))*power*power;
    double grad_psi_z = Psi_tot_*gamma*gamma*(-0.05/z*(-1.+pow(br2, power))*(-1.+pow(br1, power))*2.*power+
            0.1/z*(-1.+pow(br2, power))*(-1.+pow(br1, power))*power-
            4.*pow(br2, power)*(-1.+pow(br1, power))*power*power*
            (2.*z*Omega_*Omega_*theta_*theta_*power/(c*c*gamma*sigma_m_+z*z*Omega_*Omega_*theta_*theta_)-0.025*log(br2)/z)-
            0.1/z*pow(br1, power)*(-1.+pow(br2, power))*power*power*log(br1))/denom;
    double grad_psi_omega = Psi_tot_*8.*omega*Omega_*Omega_/c/c*pow(br1, power-1.)*(-1.+pow(br2, power))*power*power*power/denom;
    double psi = gamma*gamma*Psi_tot_*(-1.+pow(br1, power))/sigma_m_/(-1.+pow(br2, power));
    double Gamma_max = gamma + 2.*sigma_m_*psi/Psi_tot_*(1-psi/Psi_tot_);
    double epsilon = 0.5/Gamma_max/Gamma_max;
    double Omega_F = Omega_*sqrt(1.-psi/Psi_tot_);

    double E_omega = -Omega_F/2./pi/c*grad_psi_omega;
    double E_z = -Omega_F/2./pi/c*grad_psi_z;

    if (point[2] > 0) {
        return E_omega*Vector3d(cos(phi), sin(phi), 0) + E_z*Vector3d(0, 0, 1.);
    } else {
        return E_omega*Vector3d(cos(phi), sin(phi), 0) + E_z*Vector3d(0, 0, -1.);
    }
}

Vector3d BeskinMHDBField::bf_plasma_frame(const Vector3d &point, Vector3d &v) const {
   
    double x = point[0];
    double y = point[1];
    double z = point[2];
    double phi = atan2(y,x);

    double v_x = v[0];
    double v_y = v[1];
    double v_z = v[2];

    double v_omega = v_x*cos(phi)+v_y*sin(phi);
    double v_phi = v_y*cos(phi)-v_x*sin(phi);
    double v2 = v_omega*v_omega+v_phi*v_phi+v_z*v_z;

    Vector3d b = bf(point);
    double B_x = b[0];
    double B_y = b[1];
    double B_z = b[2];

    if (v2 == 0){
        return Vector3d(B_x, B_y, B_z);
    }

    if (z < 0){
        B_z = -b[2];
    }

    double B_omega = B_x*cos(phi)+B_y*sin(phi);
    double B_phi = B_y*cos(phi)-B_x*sin(phi);
    double B2 = B_omega*B_omega+B_phi*B_phi+B_z*B_z;

    Vector3d ef_ = ef(point);
    double E_x = ef_[0];
    double E_y = ef_[1];
    double E_z = ef_[2];

    double E_omega = E_x*cos(phi)+E_y*sin(phi);

    double gamma_v = 1./sqrt(1.-v2/c/c);
    double nor = (B_z*v_z+B_omega*v_omega+B_phi*v_phi)/v2;

    double B_dash_z = (gamma_v*B_z+(1-gamma_v)*v_z*nor+gamma_v/c*v_phi*E_omega)*100;
    double B_dash_omega = gamma_v*B_omega+(1-gamma_v)*v_omega*nor-gamma_v/c*v_phi*E_z;
    double B_dash_phi = gamma_v*B_phi+(1-gamma_v)*v_phi*nor+gamma_v/c*(v_omega*E_z-v_z*E_omega);

    if (z > 0) {
        return B_dash_omega*Vector3d(cos(phi), sin(phi), 0) + B_dash_phi*Vector3d(-sin(phi), cos(phi), 0) + B_dash_z*Vector3d(0, 0, 1.);
    } else {
        return B_dash_omega*Vector3d(cos(phi), sin(phi), 0) + B_dash_phi*Vector3d(-sin(phi), cos(phi), 0) + B_dash_z*Vector3d(0, 0, -1.);
    }
    
}

 /*
Vector3d BeskinMHDBField::bf_plasma_frame(const Vector3d &point, Vector3d &v) const {
    Vector3d b = bf(point);
    if (in_plasma_frame_) {
        return b;
    } else {
        return get_B_prime(b, v);
    }
}*/