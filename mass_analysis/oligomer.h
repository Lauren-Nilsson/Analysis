//
// Created by lauren on 2/22/18.
//

#ifndef LEMONSOUFFLE_OLIGOMER_H
#define LEMONSOUFFLE_OLIGOMER_H
#include <vector>
#include "vector3d.h"

class SUBUNIT;
class BEAD;

class OLIGOMER {
public:
    //data members
    unsigned int size ;
    int id ;
    double mass ;
    std::vector<SUBUNIT*> itsS;            //Gary's parent subunit
    std::vector<BEAD*> itsB;
    VECTOR3D dummy;

    OLIGOMER(VECTOR3D initial = VECTOR3D(0, 0, 0), unsigned int size_i=0, int id_i=0, double mass_i=0)       //constructor
    {
        dummy = initial;
		size = size_i;
		id = id_i;
		mass = mass_i;
    }


    ~OLIGOMER() {                                //destructor

    }

};
#endif //LEMONSOUFFLE_OLIGOMER_H
