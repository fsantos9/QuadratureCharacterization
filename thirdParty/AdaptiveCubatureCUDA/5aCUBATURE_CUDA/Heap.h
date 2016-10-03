#ifndef  HEAP_H
#define  HEAP_H

#include <algorithm>
#include <vector>
#include <ostream>

#include "CubaCuda.h"

class Heap
{
    public:
        // constructor
        Heap();

        // destructor
        //~Heap();

        void push(Region val);
        Region pop();

    private:
        std::vector<Region> values;
};




#endif  /* HEAP_H */
