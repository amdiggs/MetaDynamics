// Meta Dynamics main.cpp
// This is a mock up for the use of meta-dynamics in JDFTx to operate on a pair of ions
// No vibes, just blood,sweat, and fingertips!
// AMD 05/21/2025 

#include <cstdio>
#include <iostream>
#include <math.h>
#include "MyVec.hpp"
#include "jdftx.hpp"
#include "MetaDyn.hpp"
static Everything* e = Everything::Get();

// Computes the Periodic(closest distance between a pair of ions) in angstroms
float Periodic_Distance(AMD::Vec3 A, AMD::Vec3 B)
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
// For unit testing want to use lj force too
AMD::Vec3 Pair_LJ(AMD::Vec3 A, AMD::Vec3 B)
{
    AMD::Vec3 diff = B - A; //points from A->B
    float r = Periodic_Distance(A, B);
    //printf("r = %.3f\n",r);
    //F= -grad(U) => |F| is -1*d/dr LJ
    float mag =  12.0*(1.0/pow(r,13.0)) - 6.0*(1.0/pow(r,7.0));
    AMD::Vec3 dir = AMD::Normalize(diff);
    return mag*dir;
}

void Comp_LJ()
{
    IonInfo* inf = &(e->iInfo);
    for(int i = 0; i<inf->num_sp; i++){
	SpeciesInfo* sp = inf->m_species[i];
	for(int j = 0; j < sp->num_typ; j++){
	    AMD::Vec3 A = sp->atpos[j];
	    for(int k = j+1; k<sp->num_typ;k++){
		unsigned int a_id = i*sp->num_typ + j;
		unsigned int b_id = i*sp->num_typ + k;
		AMD::Vec3 B = sp->atpos[k];
		AMD::Vec3 F = Pair_LJ(A, B);
		//F.print();
		inf->forces[a_id] = F;
		inf->forces[b_id] = -1.0*F;
	    }
	    //additional for loop if multiple species
	}
    }
}


void Print_Pos(SpeciesInfo* spec, AMD::Vec3* orig)
{
    for(int i = 0; i<spec->num_typ; i++){
	AMD::Vec3 curr = spec->atpos[i];
	AMD::Vec3 delta = curr - orig[i];
	float d = delta.len();
	printf("Atom %d: Distance = %.3f",i,d);
    }
}

void Print_Dist(SpeciesInfo* spec)
{
    for(int i = 0; i<spec->num_typ; i++){
	AMD::Vec3 A = spec->atpos[i];
	for(int j = i+1; j<spec->num_typ; j++){
	    AMD::Vec3 B = spec->atpos[j];
	    AMD::Vec3 delta = B - A;
	    float d = delta.len();
	    printf("Atom %d-%d: Distance = %.3f\n",i,j,d);
	}
    }
}

int main(int argc, const char * argv[])
{
    e->iInfo.Init();
    SpeciesInfo* spec = e->iInfo.m_species[0];
    AMD::Vec3 init_pos[4];
    for(int i = 0; i<spec->num_typ; i++){
	init_pos[i] = spec->atpos[i];
    }
    Meta_PES mp;
    mp.Init(0.0, 7.0, 0.01, 0, 0, 0, 1);
    for(int i = 0; i<10; i++){
	Comp_LJ();
	//mp.Comp_Force();
	e->iInfo.Update_Ions();
	if(!(i%10)){
	    //Print_Pos(spec, init_pos);
	    e->iInfo.print();
	    //Print_Dist(spec);
	}
    }

    return 0;
}
