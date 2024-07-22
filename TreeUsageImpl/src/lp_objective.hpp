#ifndef LP_OJBECTIVE_HPP
#define LP_OBJECTIVE_HPP

#include <vector>
#include <string>



template<int Power>
class LPObjective
{
public:
    LPObjective(){}

    static evalutate(double value){

        return std::pow(value, Power);
    }

};

typedef LPObjective<1> KMeans;


typedef LPObjective<2> KMedian;





#endif