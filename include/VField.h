#ifndef BK_TRANSFER_VFIELD_H
#define BK_TRANSFER_VFIELD_H

#include <Eigen/Eigen>
#include "BField.h"


using Eigen::Vector3d;

// TODO: Add ``gamma`` getter - Jet instances need it if particles density in plasma frame must be calculated
class VField {
    public:
        virtual Vector3d vf(const Vector3d &point) const = 0;
        void set_profile(double r1, double r2, double K1, double K2, double b1=1.0, double b2=0.0, double b3=1.0);
    protected:
        VField(Geometry* geometry);
        // Profile parameters
        Geometry* geometry_;
        double r1_, r2_, K1_, K2_, b1_, b2_, b3_;
        bool is_profile_set;
};


class LenaMHDVField : public VField {
    public:
        LenaMHDVField(Geometry* geometry= nullptr);
        Vector3d vf(const Vector3d& point) const override;
};


class ConstFlatVField: public VField {
    public:
	    explicit ConstFlatVField(double gamma, Geometry* geometry, double betac_phi=0.0);
        Vector3d vf(const Vector3d& point) const override;

    private:
        double gamma_;
        double betac_phi_;
};


class ApproximateMHDVfield: public VField {
    public:
        ApproximateMHDVfield(Geometry* geometry, double betac_phi=0.0);
        Vector3d vf(const Vector3d& point) const override ;
    private:
        double betac_phi_;
};


class BeskinVfield: public VField {
    public:
        BeskinVfield(Geometry* geometry, BeskinMHDBField* bfield, double Psi_tot, double sigma_m, double Omega, double theta);
        Vector3d vf(const Vector3d& point) const override ;
        void set_v_phi_factor(double v_phi_factor);
    private:
        double Psi_tot_;
        double sigma_m_;
        double Omega_;
        double theta_;
        double v_phi_factor_;
        BeskinMHDBField* bfield_;
};

//
//class ShearedFlatVField: public VField {
//    public:
//        ShearedFlatVField(double gamma_axis, double gamma_border, double r);
//        Vector3d vf(const Vector3d& point) const override ;
//
//    private:
//        double gamma_axis_;
//        double gamma_border_;
//        double r_;
//};

class SheathFlatVField: public VField {
    public:
        SheathFlatVField(double gamma_spine_0, double gamma_sheath_0,
                         Geometry* geometry_in, Geometry* geometry_out,
                         double gamma_spine_v=0.0, double gamma_sheath_v=1.0,
                         double betac_phi_in=0.0,  double betac_phi_out=0.0);
        Vector3d vf(const Vector3d& point) const override ;

    private:
        double gamma_spine_0_;
        double gamma_spine_v_;
        double gamma_sheath_0_;
        double gamma_sheath_v_;
        Geometry* geometry_in_;
        double betac_phi_in_;
        double betac_phi_out_;
};

class ConstCentralVField: public VField {
    public:
	    explicit ConstCentralVField(double gamma, Geometry* geometry, double betac_phi=0.0, Vector3d origin={0, 0, 0});
        Vector3d vf(const Vector3d& point) const override;

    private:
        double gamma_;
        double betac_phi_;
        Vector3d origin_;
};

//class ShearedCentralVField: public VField {
//    public:
//        ShearedCentralVField(double gamma_axis, double gamma_border, double theta, Vector3d origin={0, 0, 0});
//        Vector3d vf(const Vector3d& point) const override ;
//
//    private:
//        double gamma_axis_;
//        double gamma_border_;
//        double theta_;
//        Vector3d origin_;
//};
//
//class SheathCentralVField: public VField {
//    public:
//        SheathCentralVField(double gamma_spine, double gamma_sheath, double theta_sheath, Vector3d origin={0, 0, 0});
//        Vector3d vf(const Vector3d& point) const override ;
//
//    private:
//        double gamma_spine_;
//        double gamma_sheath_;
//        double theta_sheath_;
//        Vector3d origin_;
//};


#endif //BK_TRANSFER_VFIELD_H
