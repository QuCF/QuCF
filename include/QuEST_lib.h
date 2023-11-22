#ifndef QUEST_LIB_H
#define QUEST_LIB_H

#include "QLib.h"


namespace YMIX{
    /** 
     * Return states with non-zero probability. 
     * @param flag_zero_ampls return zero amplitudes as well.
     * If true, bitstrings are not saved.
     * */
    void Wavefunction_Probabilities(
        const Qureg& qq, 
        StateVectorOut& out, 
        YCB flag_zero_ampls = false
    );


}

#endif



