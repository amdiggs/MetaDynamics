



#include "jdftx.hpp"
#include "MyVec.hpp"
#include <malloc/_malloc.h>
#include <random>
#include <chrono>

static unsigned int seed = (unsigned int)std::chrono::system_clock::now().time_since_epoch().count();
std::minstd_rand0 gen (seed);
// Everything/Iinfo/species
Everything* e = Everything::Get();

float urf()
{
    std::uniform_real_distribution<float> dist(0.0,1.0);
    return dist(gen);
}



SpeciesInfo::SpeciesInfo(){}

void SpeciesInfo::Init(const char* ex_symbol, float ex_rad, float ex_mass, AMD::Vec3* pos, unsigned int num)
{

    std::strcpy((char*)symbol,ex_symbol);
    rad = ex_rad;
    mass = ex_mass;
    num_typ = num;
    for (int i = 0; i<num_typ; i++){
	atpos[i] = pos[i];
    }
}


void SpeciesInfo::print()
{
    for(int i = 0; i< num_typ; i++){
	AMD::Vec3 pos = atpos[i];
	printf("%d %s: %.3f %.3f %.3f\n",i, symbol, pos.x,pos.y,pos.z);
    }
}



IonInfo::IonInfo():num_sp(0), init(false){}

void IonInfo::Init()
{
    this->m_coordsType = e->g_coords;
    m_species = (SpeciesInfo**)malloc(sizeof(SpeciesInfo*));
    m_species[0] = (SpeciesInfo*)malloc(sizeof(SpeciesInfo));
    AMD::Vec3 tmp[4] = {AMD::Vec3(0.0,0.0,0.0), AMD::Vec3(0.4,0.0,0.0),AMD::Vec3(0.0,0.8,0.0), AMD::Vec3(-0.8,-0.8,0.0)};
    m_species[0]->Init("O", 0.65,15.998,tmp,4);
    num_sp = 1;

    for(int i = 0; i < num_sp; i++){
	SpeciesInfo* spec = m_species[i];
	for(int j = 0; j < spec->num_typ; j++){
	    unsigned int f_id = i*spec->num_typ + j;
	    forces[f_id] = 0.;
	    velocities[f_id] = 0.; //0.1f*urf();
	}
    }
    init = true;
}

IonInfo::~IonInfo()
{
    if(init){
	for(int i =0; i<num_sp; i++){
	    free(m_species[i]);
	}
	free(m_species);
    }
}


unsigned int IonInfo::Get_SP_Indx(const char* type)
{
    return 0;
}


void IonInfo::print(){
    for(int i = 0; i < 4; i++){
	AMD::Vec3 f = forces[i]; 
	printf("Force on Atom %d: (%.3f, %.3f,%.3f)\n",i,f.x,f.y,f.z);
    }
}


void IonInfo::Update_Ions()
{
    float dt = e->dt;
    for(int i = 0; i < num_sp; i++){
	SpeciesInfo* spec = m_species[i];
	for(int j = 0; j < spec->num_typ; j++){
	    unsigned int f_id = i*spec->num_typ + j;
	    AMD::Vec3 a = 1.0*forces[f_id];
	    //a.print();
	    velocities[f_id] += dt*a;
	    AMD::Vec3 dx = dt*velocities[f_id];
	    printf("######Atom %d########\n",j);
	    spec->atpos[j].print();
	    spec->atpos[j] += dx;
	    dx.print();
	    spec->atpos[j].print();
	}
    }
    e->time+=dt;
}


Everything::Everything()
{
    lattice[0] = AMD::Vec3(6.5,0.0,0.0);
    lattice[1] = AMD::Vec3(0.0,6.5,0.0);
    lattice[2] = AMD::Vec3(0.0,0.0,21.2);
    LattM.assign_row(0, lattice[0]);
    LattM.assign_row(1, lattice[1]);
    LattM.assign_row(2, lattice[2]);
    invR = AMD::Inverse(LattM);
    g_coords = CoordsCartesian;
}

Everything Everything::inst;

Everything* Everything::Get(){return &inst;}


