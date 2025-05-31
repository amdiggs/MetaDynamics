//
//
//
#ifndef JDFTX_H
#define JDFTX_H


#include "stdio.h"
#include "MyVec.hpp"

//! Coordinate system for ion positions
enum CoordsType {CoordsLattice, CoordsCartesian}; 

//! Coordinate system for force output:
enum ForcesOutputCoords { ForcesCoordsPositions, ForcesCoordsLattice, ForcesCoordsCartesian, ForcesCoordsContravariant };

struct AtomicMode
{	unsigned int sp, at; //!< species and atom number
	unsigned int Force_id;
	AMD::Vec3 dirLattice; //!< Perturbation in lattice coords
	AMD::Vec3 dirCartesian; //!< Perturbation in Cartesian coords
};

struct SpeciesInfo
{
    char symbol[2];
    float rad;
    float mass;
    AMD::Vec3 atpos[10];
    unsigned int num_typ;
    SpeciesInfo();
    void Init(const char* symbol, float rad, float mass, AMD::Vec3* pos, unsigned int num_ats);
    ~SpeciesInfo(){};
    void print();
};



struct IonInfo
{

    unsigned int num_sp;
    SpeciesInfo** m_species = NULL;
    CoordsType m_coordsType;
    AMD::Vec3 forces[10];
    AMD::Vec3 velocities[10];
    IonInfo();
    IonInfo(const IonInfo&) = delete;
    IonInfo& operator=(const IonInfo&) = delete;
    void Init();
    ~IonInfo();
    bool init;
    unsigned int Get_SP_Indx(const char*);
    void print();
    void Update_Ions();
};

class Everything
{
private:
    AMD::Vec3 lattice[3];
    Everything();
    static Everything inst;
public:
    Everything(Everything&) = delete;
    void operator=(const Everything&) = delete;
    static Everything* Get();
    AMD::Mat3 LattM;
    AMD::Mat3 invR;
    IonInfo iInfo;
    CoordsType g_coords;
    float time = 0.;
    float dt = 0.1;

};

#endif
