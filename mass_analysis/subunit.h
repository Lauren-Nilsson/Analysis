//
// Created by lauren on 1/25/18.
//

#ifndef LEMONSOUFFLE_UNIT_H
#define LEMONSOUFFLE_UNIT_H

#include <vector>
#include "vector3d.h"
#include "oligomer.h"


class BEAD;


class SUBUNIT
{
public:
    //member variables
    int id;                                 //subunit id
    std::vector<BEAD*> itsB;                //particles making up the subunit
    std::vector<OLIGOMER> itsO;                         //its oligmer
   
    VECTOR3D vsumvec;


    //member fxns
    SUBUNIT(VECTOR3D initial=VECTOR3D(0,0,0))       //constructor
    {
        vsumvec=initial;
    }




};




#endif //LEMONSOUFFLE_UNIT_H
