

#include "MetaDyn.hpp"
#include "jdftx.hpp"


extern Everything* e;



static float Periodic_Distance(AMD::Vec3 A, AMD::Vec3 B)
{
    CoordsType c = e->g_coords;
    if(c == CoordsCartesian){
	A = e->invR*A;
	B = e->invR*B;
    }
    AMD::Vec3 diff = B-A;

    for(int i = 0; i<3;i++){
	float x = (abs(diff[i]) > 0.5) ? 1.0 - abs(diff[i]) : abs(diff[i]);
	diff[i] = x;
    }
    diff = e->LattM*diff;
    return diff.len();
}

//This is a cheat of a gaussian where we will only use delta x because mu is the current value of the CV
float Gauss_of_Energy(float del_x, float sigma)
{
    float N = 1.0 / sqrt(twoPI*sigma*sigma);
    float arg = (del_x*del_x)/(2.0*sigma*sigma);
    return N*exp(-1.0*arg);
}


Meta_PES::Meta_PES(){}
Meta_PES::~Meta_PES()
{
    if(init){
	free(m_potential);
    }
}

void Meta_PES::Init(float  low, float high, float dx, unsigned int sp_A, unsigned int at_A, unsigned int sp_B, unsigned int at_B)
{
    this->m_low = low;
    this->m_high = high;
    this->m_dx = dx;
    this->m_sigma = dx;
    this->ion_A.sp = sp_A;
    this->ion_A.at = at_A;
    this->ion_B.sp = sp_B;
    this->ion_B.at = at_B;
    unsigned int A_idx = 0;
    unsigned int B_idx = 0;
    //I added a member var to Atomic mode which for quick access to the force index in 
    //IonInfo, There is probably a much better way!!
    IonInfo* local_Ions = &(e->iInfo);
    for(int i = 0;i< local_Ions->num_sp;i++){
	SpeciesInfo* spec = local_Ions->m_species[i];
	if(i<ion_A.sp){A_idx+=spec->num_typ;}
	if(i<ion_B.sp){B_idx+=spec->num_typ;}
    }
    ion_A.Force_id = A_idx + ion_A.at;
    ion_B.Force_id = B_idx + ion_B.at;

    m_num_bins = floor((high - low)/dx);
    m_potential = (float*)malloc(4*m_num_bins);
    for (int i =0; i< m_num_bins; i++){
	m_potential[i] = 0.;
    }
    init = true;
}


void Meta_PES::Comp_Force()
{
    SpeciesInfo* local_SP_A = e->iInfo.m_species[ion_A.sp];
    SpeciesInfo* local_SP_B = e->iInfo.m_species[ion_B.sp];
    AMD::Vec3 A_pos = local_SP_A->atpos[ion_A.at];
    AMD::Vec3 B_pos = local_SP_B->atpos[ion_B.at];
    float dist = Periodic_Distance(A_pos, B_pos); //this distance is gurenteeded to be in A
    int bin = floor(dist / m_dx);
    if(bin > (m_num_bins-3)){return;}
    int low = (bin >= 2) ? -2 : -bin;

    // here we will add a Gaussian dolup of energy centered at the current reaction coordinate.
    for(int i = low; i<3; i++){
	int idx = bin + i;
	m_potential[idx] += Gauss_of_Energy((float)idx*m_dx, m_sigma);
    }
    float der =  m_potential[bin + 1] - m_potential[bin - 1] ;
    der /=(2.0*m_dx);
    // A positive F along a reaction coordinate is in our case pusing the atoms further apart
    AMD::Vec3 dir = AMD::Normalize(B_pos - A_pos);
    //VB - VA points from A -> B
    // I don't think there need to be a coordstype check but I could be wrong
    float scale = 1.0; // this should be some appropriate energy scale based on how many steps we want
    AMD::Vec3 F_A = -1.0*der*scale*dir;
    AMD::Vec3 F_B = -1.0*F_A;
    e->iInfo.forces[ion_A.Force_id] += F_A;
    e->iInfo.forces[ion_B.Force_id] += F_B;
return;
}


