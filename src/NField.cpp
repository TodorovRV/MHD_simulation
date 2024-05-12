#include "utils.h"
#include "NField.h"
#include "MyExceptions.h"
#include "BField.h"


NField::NField(bool in_plasma_frame, ParticlesDistribution* particles, Geometry* geometry_out, Geometry* geometry_in,
			   VField* vfield) :
    in_plasma_frame_(in_plasma_frame) {
    particles_ = particles;
    geometry_in_ = geometry_in;
    geometry_out_ = geometry_out;
	vfield_ = vfield;
}

double NField::nf(const Vector3d &point) const {
    double x, y, r_point, r_border, result;
    x = point[0];
    y = point[1];
    r_point = sqrt(x*x + y*y);
	
	
	double factor = 1.0;
	
	if(geometry_out_) {
		r_border = geometry_out_->radius_at_given_distance(point);
		if (r_point > r_border) {
			factor = exp(-pow(r_point - r_border, 2.0)/l_eps_N/l_eps_N);
		}
	}
//	if(geometry_in_) {
//		r_border = geometry_in_->radius_at_given_distance(point);
//		if (r_point < r_border) {
//			factor = exp(-pow(r_point - r_border, 2.0)/l_eps_N/l_eps_N);
//		}
//	}
	return _nf(point)*factor;
	
//    if(geometry_out_) {
////        std::cout << "Zeroing N' because of the OUTER geometry!" << "\n";
//        // Find radius of outer border at given point
//        r_border = geometry_out_->radius_at_given_distance(point);
//        if (r_point > r_border) {
//            return 0.0;
//        }
//    }
//    if(geometry_in_) {
////        std::cout << "Zeroing N' because of the INNER geometry!" << "\n";
//        // Find radius of inner border at given point
//        r_border = geometry_in_->radius_at_given_distance(point);
//        if (r_point < r_border) {
//            return 0.0;
//        }
//    }
//    return _nf(point);
}

double NField::nf_plasma_frame(const Vector3d &point, double &gamma) const {
    double n = nf(point);
    if (in_plasma_frame_) {
        return n;
    } else {
        return n/gamma;
    }
}


BKNField::BKNField(double n_0, double n_n, ParticlesDistribution* particles, bool in_plasma_frame,
                   Geometry* geometry_out, Geometry* geometry_in,  VField* vfield) :
        NField(in_plasma_frame, particles, geometry_out, geometry_in, vfield),
        n_0_(n_0),
        n_n_(n_n),
        r_mean_(0.0),
        r_width_(0.0),
        background_fraction_(0.0),
        is_profile_set_(false) {}


double BKNField::_nf(const Vector3d &point) const {
    double r = point.norm();
    double raw_density = n_0_ * pow(r / pc, -n_n_);
    if(is_profile_set_){
        double x = point[0];
        double y = point[1];
        double R_cur = hypot(x, y);
        double R_out = geometry_out_->radius_at_given_distance(point);
//        if(R_cur/R_out < 0.25){
//            std::cout << R_cur/R_out << "\n";
//        }
        return (generalized1_gaussian1d(R_cur / R_out, r_mean_, r_width_, 10) + background_fraction_) * raw_density;
    }else{
        return raw_density;
    }
}


void BKNField::set_heating_profile(double r_mean, double r_width, double background_fraction) {
    if (r_mean > 1 || r_mean < 0.0 || r_width < 0.0 || background_fraction < 0.0 || background_fraction_ > 1.0){
        throw BadHeatingProfileParameters();
    }
    r_mean_ = r_mean;
    r_width_ = r_width;
    background_fraction_ = background_fraction;
    is_profile_set_ = true;
}


GJNfield::GJNfield(double n_0, double Psi_tot, double sigma_m, double Omega, double theta,
                ParticlesDistribution* particles, bool in_plasma_frame, Geometry* geometry_out, Geometry* geometry_in,  VField* vfield) :
    NField(in_plasma_frame, particles, geometry_out, geometry_in, vfield),
    n_0_(n_0), 
    Psi_tot_(Psi_tot), 
    sigma_m_(sigma_m), 
    Omega_(Omega), 
    theta_(theta),
    r_mean_(0.0),
    r_width_(0.0),
    background_fraction_(0.0),
    is_profile_set_(false)  {}

double GJNfield::_nf(const Vector3d &point) const {
    double x = point[0];
    double y = point[1];
    double z = point[2];
    double phi = atan2(y,x);
    double gamma = 1;
    double omega = sqrt(x*x+y*y);
    double power = 1-0.025*log(z*abs(Omega_)*theta_/c/sqrt(gamma*sigma_m_));
    double br1 = 1+(omega*omega*Omega_*Omega_/c/c/gamma/gamma);
    double br2 = 1+(z*z*Omega_*Omega_*theta_*theta_/c/c/gamma/sigma_m_);
    double denom = 4*sigma_m_*(-1+pow(br2, power))*(-1+pow(br2, power))*power*power;
    double psi = gamma*gamma*Psi_tot_*(-1+pow(br1, power))/sigma_m_/(-1+pow(br2, power));

    double grad_psi_omega = Psi_tot_*8*omega*Omega_*Omega_/c/c*pow(br1, power-1)*(-1+pow(br2, power))*power*power*power/denom;
    double Omega_F = Omega_*sqrt(1-psi/Psi_tot_);

    double B_z = grad_psi_omega/2/pi/omega;

    double raw_density = n_0_*abs(B_z*Omega_F)/2/pi/c/q_e;

    if(is_profile_set_){
        double x = point[0];
        double y = point[1];
        double R_cur = hypot(x, y);
        double R_out = geometry_out_->radius_at_given_distance(point);
        return (generalized1_gaussian1d(R_cur / R_out, r_mean_, r_width_, 10) + background_fraction_) * raw_density;
    } else {
        return raw_density;
    }
}

void GJNfield::set_heating_profile(double r_mean, double r_width, double background_fraction) {
    if (r_mean > 1 || r_mean < 0.0 || r_width < 0.0 || background_fraction < 0.0 || background_fraction_ > 1.0){
        throw BadHeatingProfileParameters();
    }
    r_mean_ = r_mean;
    r_width_ = r_width;
    background_fraction_ = background_fraction;
    is_profile_set_ = true;
}


B2Nfield::B2Nfield(double n_0, double Psi_tot, double sigma_m, double Omega, double theta, std::vector<VectorBField*> vbfields,
                ParticlesDistribution* particles, bool in_plasma_frame, Geometry* geometry_out, Geometry* geometry_in,  VField* vfield) :
    NField(in_plasma_frame, particles, geometry_out, geometry_in, vfield),
    n_0_(n_0), 
    Psi_tot_(Psi_tot), 
    sigma_m_(sigma_m), 
    Omega_(Omega), 
    theta_(theta),         
    r_mean_(0.0),
    r_width_(0.0),
    background_fraction_(0.0),
    is_profile_set_(false) 
    {vbfields_ = vbfields;}

double B2Nfield::_nf(const Vector3d &point) const {
    // TODO: this is ok only for one Bfield
    auto v = vfield_->vf(point);
    double b = 0.0;
    for(auto bfield: vbfields_) {
		Vector3d vb = bfield->bf_plasma_frame(point, v);
		b += vb.norm();
	}
    double raw_density = n_0_*b*b;
    if(is_profile_set_){
        double x = point[0];
        double y = point[1];
        double R_cur = hypot(x, y);
        double R_out = geometry_out_->radius_at_given_distance(point);
//        if(R_cur/R_out < 0.25){
//            std::cout << R_cur/R_out << "\n";
//        }
        return (generalized1_gaussian1d(R_cur / R_out, r_mean_, r_width_, 10) + background_fraction_) * raw_density;
    }else{
        return raw_density;
    }
}

void B2Nfield::set_heating_profile(double r_mean, double r_width, double background_fraction) {
    if (r_mean > 1 || r_mean < 0.0 || r_width < 0.0 || background_fraction < 0.0 || background_fraction_ > 1.0){
        throw BadHeatingProfileParameters();
    }
    r_mean_ = r_mean;
    r_width_ = r_width;
    background_fraction_ = background_fraction;
    is_profile_set_ = true;
}


j2Nfield::j2Nfield(double n_0, double Psi_tot, double sigma_m, double Omega, double theta, std::vector<VectorBField*> vbfields,
                ParticlesDistribution* particles, bool in_plasma_frame, Geometry* geometry_out, Geometry* geometry_in,  VField* vfield) :
    NField(in_plasma_frame, particles, geometry_out, geometry_in, vfield),
    n_0_(n_0), Psi_tot_(Psi_tot), sigma_m_(sigma_m), Omega_(Omega), theta_(theta) {vbfields_ = vbfields;}

double j2Nfield::_nf(const Vector3d &point) const {
    double x = point[0];
    double y = point[1];
    double z = point[2];
    double phi = atan2(y,x);
    double omega = sqrt(x*x+y*y);
    double d = 0.0001;

    auto v = vfield_->vf(point);
    double j2 = 0;

    // TODO: this is ok only for one Bfield
    // Due to comlexity of equations approximating derivative as a differential difference
    // Using notion: B<1 or 0>_<direction of derivative>_<comonent>
    // 0 for value in the <point>, 1 for value in <point> + d<point>
    for(auto bfield_: vbfields_) {
        auto b0 = bfield_->bf_plasma_frame(point, v);
        double B0_x = b0[0];
        double B0_y = b0[1];
        double B0_z = b0[2];
        double B0_omega = B0_x*cos(phi)+B0_y*sin(phi);
        double B0_phi = B0_y*cos(phi)-B0_x*sin(phi);

        auto b1_z = bfield_->bf_plasma_frame(Vector3d(omega*cos(phi), omega*sin(phi), z+d*z), v);
        double B1_z_x = b1_z[0];
        double B1_z_y = b1_z[1];
        double B1_z_z = b1_z[2];
        double B1_z_omega = B1_z_x*cos(phi)+B1_z_y*sin(phi);
        double B1_z_phi = B1_z_y*cos(phi)-B1_z_x*sin(phi);    

        auto b1_omega = bfield_->bf_plasma_frame(Vector3d((omega+d*omega)*cos(phi), (omega+d*omega)*sin(phi), z), v);
        double B1_omega_x = b1_omega[0];
        double B1_omega_y = b1_omega[1];
        double B1_omega_z = b1_omega[2];
        double B1_omega_omega = B1_omega_x*cos(phi)+B1_omega_y*sin(phi);
        double B1_omega_phi = B1_omega_y*cos(phi)-B1_omega_x*sin(phi);   

        auto b1_phi = bfield_->bf_plasma_frame(Vector3d(omega*cos(phi+d*phi), omega*sin(phi+d*phi), z), v);
        double B1_phi_x = b1_phi[0];
        double B1_phi_y = b1_phi[1];
        double B1_phi_z = b1_phi[2];
        double B1_phi_omega = B1_phi_x*cos(phi)+B1_phi_y*sin(phi);
        double B1_phi_phi = B1_phi_y*cos(phi)-B1_phi_x*sin(phi);    

        double j_z = c/4./pi/omega*(((omega+d*omega)*B1_omega_phi - omega*B0_omega)/(d*omega) - (B1_phi_omega - B0_omega)/(d*phi));
        double j_omega = c/4./pi*((B1_phi_z - B0_z)/(d*phi)/omega - (B1_z_phi - B0_phi)/(d*z));
        double j_phi = c/4./pi*((B1_z_omega - B0_omega)/(d*z) - (B1_omega_z - B0_z)/(d*omega));
        j2 += j_z*j_z+j_omega*j_omega+j_phi*j_phi;
    }
    return n_0_*j2;
}


EquipartitionBKNfield::EquipartitionBKNfield(ParticlesDistribution *particles,
											 std::vector<VectorBField*> vbfields,
											 Geometry *geometry_out, Geometry *geometry_in, VField *vfield,
											 double fraction):
	NField(true, particles, geometry_out, geometry_in, vfield),
	fraction_(fraction) {
	vbfields_ = vbfields;
}

double EquipartitionBKNfield::_nf(const Vector3d &point) const {
	auto v = vfield_->vf(point);
	double b = 0.0;
	for(auto bfield: vbfields_) {
		Vector3d vb = bfield->bf_plasma_frame(point, v);
		b += vb.norm();
	}
    bool set_sheath = false;
    double raw_density = fraction_*particles_->get_equipartition_bsq_coefficient()*b*b;
    double background_fraction = 1.0;
    double k_ = 1.0;
    double amp = 0.005*pc;
    double width = 0.05*pc;
    double impact_factor = 10;
    double density = 0;
    double x = point[0];
    double y = point[1];
    double z = point[2];
    if (set_sheath){
        double distance = amp*pow(abs(z)/pc, k_);
        density = raw_density*(background_fraction+impact_factor*
                generalized1_gaussian1d(hypot(x, y)/pc, distance/pc, width/pc*pow(abs(z)/pc, k_), 2));
    } else {
        density = raw_density;
    }
	return density;
}
