// This is a meta-dynamics class for use in JDFTx
// Andrew Diggs 05-31-2025
#ifndef META_DYNAMICS
#define META_DYNAMICS


#include "stdio.h"
#include "jdftx.hpp"

// For the meta dynamics we want to define a collective variable (CV).
// The CV will be asscociated with a potential energy surface (PES).
// For our case of O-O bond breaking we will use the Bond length.
// Thus the surface will be a curve => PES -> PEC
//

class Meta_PES
{
private:
    // float m_low is the lower cuttoff for the PES. We can set m_low to 0 for ease of use,
    // the atoms should never be 0 angstrom apart so to optimize.
    // float high will be the upper cutoff, at ~5-6 A the O-O bond is throughly broken
    // We can use 6 to be safe but better optimization would reduce memory cost.
    float m_low, m_high, m_dx, m_sigma;
    //CV_func m_CV;
    // num bins will be the numer of descrete values repersenting our PES, this will be given by (high - low) / dx
    unsigned int m_num_bins;
    float* m_potential = NULL;
    bool init = false;
    //AtomicMode was taken from JDFTX>Perturb
    //struct that contains int species num, int ion index vec3 F in cartesin and vec3 F in lattice coords 
    // to be used with Ioninfo
    AtomicMode ion_A, ion_B;
    // Because of the dynamic allocation of m_potential I am just going to delete the copy and clone constructors
    // If a deep copy is better I can redo this
    Meta_PES(const Meta_PES&) = delete;
    Meta_PES operator=(const Meta_PES&) = delete;
public:
    Meta_PES();
    ~Meta_PES();
    void Init(float low, float high, float dx, unsigned int sp_A, unsigned int at_A, unsigned int sp_B, unsigned int at_B);
    void Comp_Force();
};

#endif // !META_DYNAMICS
