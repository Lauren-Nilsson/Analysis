//
// Created by lauren on 1/25/18.
//

#ifndef LEMONSOUFFLE_BEAD_H
#define LEMONSOUFFLE_BEAD_H

#include <vector>
#include <cmath>
#include "vector3d.h"



class SUBUNIT;
class OLIGOMER;



class BEAD
{
public:

//member variables
    int id ;                         //particle ID
    double m;                           //mass of the particle (amu)
    double sigma;                       //diameter of particle (unitless)
    VECTOR3D bx;
    int unit ;
    int type ;
    VECTOR3D pos;                       //position of the particle (xyz) (unitless)
    VECTOR3D vel;                       //velocity of the particle (xyz) (delt^-1)
    std::vector<SUBUNIT*> itsS;        
    std::vector<OLIGOMER> itsSO;            //its sub_oligomer
    OLIGOMER* itsO;
	
 


//member functions

    BEAD(VECTOR3D position_i=VECTOR3D(0,0,0),double mi=0, int id_i=0, int unit_i=0, int type_i=0, OLIGOMER* itsO_i=NULL)
    {                                                       //constructor fxn
        m = mi;
        pos = position_i;
		id = id_i;
		unit = unit_i;
		type = type_i;
		itsO = itsO_i;
    }

  


};


#endif //LEMONSOUFFLE_BEAD_H
