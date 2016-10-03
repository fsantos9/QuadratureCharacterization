#include "Heap.h"


// constructor
Heap::Heap()
{ }



static bool cmpMax(const Region &a, const Region &b) {
    return a.ErrorMax < b.ErrorMax;
}


void Heap::push(Region val)
{
    values.push_back(val);
    push_heap(values.begin(), values.end(), cmpMax);
}


Region Heap::pop()
{
    pop_heap(values.begin(), values.end(), cmpMax);
    Region val = values.back();
    values.pop_back();
    return val;
}
