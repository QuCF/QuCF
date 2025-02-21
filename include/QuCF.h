#ifndef ORACLETOOL_H
#define ORACLETOOL_H

#include "BaseTool.h"

class QuCF__ : public BaseTool__
{
public:
    /**
     * @param[in] project_name  is a project name that defines names of input and output files.
     * @param[in] path_to_inputs is path to input files.
     */
    QuCF__(
        const QuESTEnv& env, 
        YCS project_name, 
        YCS path_to_inputs 
    );
    ~QuCF__();
    void launch();

protected:
    void read_circuit_structure_from_file(YCS data);
    

private:
    void read_constants(YISS str);
    void read_options(YISS str);
    void read_circuit_declaration(YISS istr);
    void read_circuit_structure(YISS istr, YSQ* oc_ext = nullptr);
    void read_main_circuit(YISS istr);
    void read_gate(YISS istr, YPQC oc, YCB flag_inv=false);

    /**
     * @param oc: parent circuit;
     * @param flag_inv: if true, inverse subcircuit will be inserted to the parent circuit.
    */
    void read_subcircuit(YISS istr, YPQC oc, YCB flag_inv=false);

    /**
     * Read gates from another .oracle file
    */
    void read_gates_from_file(YISS istr, YPQC oc);

    void read_state(YISS istr);
    void read_state_init_file();
    qreal get_value_from_word(YCS word);
    void  calc(std::shared_ptr<QCircuit>& u_work, YCI count_init_state);

    /**
     * sel_matrix = 1: the matrix entangled with the zero state of all ancillary qubits;
     * sel_matrix = 1: the entire unitary matrix of the circuit.
     * Here, it is assumed that all ancillary qubits have a higher significance than the state qubits.
    */
    void calc_matrix(
        std::shared_ptr<QCircuit>& u_work, 
        bool flag_sparse_format=false,
        int sel_matrix=1
    );

private:
    // dictionary with constants to create the oracle:
    std::map<std::string, qreal> constants_; 

    // several initial states:
    // every init. state is represented by several registers, 
    //  where several qubits might be set to 1.
    std::vector<std::map<std::string, std::vector<int>>> init_states_;

    // amplitudes of the initial state read from the initial state file:
    std::vector<qreal> init_ampl_vec_real_;
    std::vector<qreal> init_ampl_vec_imag_;

    // if true, the initial state is read from the .init_state file:
    bool flag_init_state_file_;

    /**
     * true: construct a matrix, where
     * positions of matrix elements are indicated by nonancilla state bitstrings,
     * and elements' values are indicated by the state amplitudes
     * (only states entangled with zero-ancilla states are considered);
     * the resulting matrix is written to the .hdf5 file;
     */ 
    int flag_matrix_;

    /**
     * none          - do not compute/print any output states (no printing on the sceen);
     * all           - compute all output states (and separately, zero-ancilla states);
     *                     however, do not print separately the zero-ancilla states;
     * zero-ancillae - compute/print states where all ancillae are in the zero state
     * (remark, ancilla qubits are assumed to have higher priority than the main qubits).
     */ 
    std::string sel_compute_output_; 
    std::string sel_print_output_; 

    /**
     * To calculate probabilities of all states for the given qubits, \p focus_qubits_.
     */
    bool flag_prob_;
    YVIv focus_qubits_;

    /**
     * Output data from QSP
    */
    bool flag_output_qsp_;

    /**
     * Output data from gadget.
    */
    bool flag_output_gadget_;

    /**
     * True: include stop gates into circuits.
    */
    bool flag_stop_gates_;

    /**
     * True: show a progress bar if flag_stop_gates_ = False
     */
    bool flag_progress_bar_;

    /**
     * True: in the repeat gate, simply insert gates.
    */
    bool flag_repeat_insert_;

    /**
     * Output amplitudes of all states entagled with the zero state of the following qubits.
    */
    std::vector<short>  qs_chosen_state_;
    bool flag_zero_state_of_;

    // Information about complex structures within the circuit:
    QuCF_complex_data cd_;
};
#endif

