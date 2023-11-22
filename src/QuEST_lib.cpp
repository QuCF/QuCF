#include "../include/QuEST_lib.h"

using namespace std;


void YMIX::Wavefunction_Probabilities(
    const Qureg& qq, 
    StateVectorOut& out, 
    YCB flag_get_zero_ampls
){ 
    const unsigned n = qq.numQubitsRepresented; // number of qubits;
    long long N = pow(2, out.n_low_prior_qubits);
    Complex aa;
    qreal prob;
    bool flag_chosen;
    bool flag_to_choose; 

    // check whether it is necessary to choose a special state or not
    flag_to_choose = true;
    if(out.state_to_choose.empty()) flag_to_choose = false;

    // find states available for such a number of qubits and their amplitudes:
    out.ampls.clear();
    out.states.clear();
    copyStateFromGPU(qq);
    for(unsigned id_state = 0; id_state < N; id_state++)
    {
        // aa = getAmp(qq, id_state);
        aa.real = qq.stateVec.real[id_state];
        aa.imag = qq.stateVec.imag[id_state];

        // check whether the corresponding probability is not equal to zero
        prob = sqrt(pow(aa.real, 2) + pow(aa.imag, 2));
        if(!YMATH::is_zero(prob) or flag_get_zero_ampls)
        {
            flag_chosen = true;
            vector<short> one_state(n);
            YMATH::intToBinary(id_state, one_state);
            if(flag_to_choose)
                for(unsigned id_qubit = 0; id_qubit < out.state_to_choose.size(); ++id_qubit)
                {
                    if(out.state_to_choose[id_qubit] < 0)
                        continue;
                    if(out.state_to_choose[id_qubit] != one_state[id_qubit])
                    {
                        flag_chosen = false;
                        break;
                    }
                }
            if(flag_chosen)
            {
                out.ampls.push_back(aa);
                if(not flag_get_zero_ampls)
                    out.states.push_back(one_state);
            }
        }
    }
    getStrWavefunction(out);
}


