#ifndef LP_OJBECTIVE_HPP
#define LP_OBJECTIVE_HPP

#include <cmath>


template<int Power>
class LPObjective
{
public:
    LPObjective(){}

    static double Evaluate(double value){

        return std::pow(value, Power);
    }

};

typedef LPObjective<1> KMedian;


typedef LPObjective<2> KMeans;





#endif