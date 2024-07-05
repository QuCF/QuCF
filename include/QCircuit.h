#ifndef QCIRCUIT_H
#define QCIRCUIT_H


#include "QGates.h"
#include "CircuitLayers.h"

/**
 * @brief Circuit class.
 * ---
 * The 0-th qubit is the least significant 
 * (the rightmost in the statevector and the very bottom qubit in the quantum circuit)
 */
class QCircuit{
    public:
    /** Constructor of a circuit
     * @param[in] name name of the circuit;
     * @param[in] env input QuEST execution environment; 
     * @param[in] path_to_output path where output files should be written;
     * @param[in] nq number of qubits to create. 
     * If \p nq <= 0, the QuEST circuit will not be created. 
     * In this case, it will be necessary to launch separately the function "create";
     * @param[in] constants dictionary of constans used to create the circuit (by default, empty);
     * @param[in] flag_circuit print or not .circuit files (by default, false);
     * @param[in] flag_tex print or not .tex files (by default, false);
     * @param[in] flag_layers to calculate or not the layer for each gate;
     * */
    QCircuit(
        YCS name, const QuESTEnv& env, YCS path_to_output="./", YCU nq = 0,
        const std::map<std::string, qreal>& constants = std::map<std::string, qreal>(),
        YCB flag_circuit = false,
        YCB flag_tex = false,
        YCB flag_layers = false,
        YCB flag_stop_gates = true,
        YCB flag_repeat_insert = false
    );

    /**
     * @brief Copy a circuit.
     * @param[in] oc is a circuit to copy.
     * @param[in] cname new name of the circuit.
     */
    QCircuit(YCCQ oc, YCS cname="");

    /** Destructor */
    ~QCircuit();

    /**
     * @brief Create a QuEST circuit.
     * @param nq number of qubits to create (has to be > 0);
     */
    void create(YCU nq);
    void allocate_circuit();

    void create_circ_file();
    void create_tex_file();
    void finish_tex_file();

    void print_gates();

    inline std::string get_name(){ return name_; }
    inline unsigned get_n_qubits() const { return nq_; }
    inline int get_n_gates(){ return gates_.size(); }
    std::map<std::string, YVIv> get_regs(){ return regs_; }
    YVSv get_reg_names() const { return regnames_; }
    unsigned get_nq_in_reg(YCS reg_name) const{ return regs_.at(reg_name).size(); }

    /**
     * @brief Generate the circuit.
     */
    void generate();

    /**
     * @brief Generate the circuit taking into account stop gates.
     */
    void generate(std::string& stop_name, int& id_current);

    void activate_gadget(YCI id_counter, YCI N_mult)
    {
        auto start = gates_.begin();
        for(auto it = start; it != gates_.end(); ++it)
            (*it)->activate_gadget(id_counter, N_mult);
    }
    void deactivate_gadget();

    /**
     * @brief Transform the current matrix into its Hermitian adjoint version.
     */
    void h_adjoint();


    /**
     * @brief Add control qubits to the circuit.
     * @param[in] cs_unit 1-control nodes to add.
     * @param[in] cs_zero 0-control nodes to add.
     */
    inline void add_control_qubits(YCVI cs_unit, YCVI cs_zero = {})
    {
        for(auto& gate: c->gates_)
        {
            gate->add_control_qubits(cs_unit, cs_zero);
        }
    }



    inline Qureg get_qureg(){return c_;}

    /**
     * @brief Copy gates from the circuit \p circ, correct their positions \p regs_new, 
     * and add them to the end of the current circuit.
     * @param circ circuit from where new gates are taken;
     * @param regs_new registers of the current circuit to which the new gates are applied.
     *  E.g., \p circ has two qubits: 0, 1 (0 - least significant (bottom) qubit).
     *  If \p regs_new = [3,2], then 0->3, 1->2: zero qubit of \p circ is applied to 
     *  the qubit 3 of the current circuit, while qubit 1 of \p circ is applied to
     *  the qubit 2 of the current circuit.
     * @param box a ghost gate to indicate boundaries where 
     *  gates from the circuit \p circ are placed into the current circuit.
     * @param flag_inv if true, first get adjoint gates from \p circ. 
     * @param cs_unit unit-control qubits that should control each gate from \p circ.
     * @param cs_zero zero-control qubits that should control each gate from \p circ.
     */
    void copy_gates_from(
        YCCQ circ, 
        YCVI regs_new, 
        YCCB box = std::shared_ptr<const Box__>(nullptr),
        YCVI cs_unit = YVIv {},
        YCVI cs_zero = YVIv {},
        YCB flag_inv = false
    );

    /**
     * @brief Insert gates (!!! without copying them and without the correction of their qubits !!!) 
     * from the circuit \p of to the end of the current circuit.
     * @param of circuit from where new gates are taken;
     */
    void insert_gates_from(
        const QCircuit* of, 
        YCCB box = std::shared_ptr<const Box__>(nullptr)
    );

    /**
     * @brief Add a register. 
     * The first added register is placed at the top of the circuit 
     * (is the most significant register).
     * @param[in] name name of the register;
     * @param[in] n_qubits number of qubits in the register;
     * @param[in] flag_ancilla if true, it is an ancilla register;
     * @return qubits of the register. 
     * The 0-th element in the returned vector corresponds to the least significant qubit.
     */
    YVIv add_register(YCS name, YCU n_qubits, YCB flag_ancilla=false);

    void set_standart_output_format();

    inline void set_ancillae(YCVI ancs){ ancs_ = YVIv(ancs); }

    /**
     * @brief Get the ancillae qubits of the circuit.
     * The vector with ancillae is not ordered.
     * @return YVIv unordered vector with ancillae qubits of the circuit.
     */
    inline YVIv get_ancillae() const { return ancs_; }

    /**
     * @brief Get an output format that corresponds to 
     *        the sizes of implied registers of the circut.
     * @return YVUv - vector that represents the format.
     */
    inline YVUv get_standart_output_format(){ return standart_output_format_; }

    /**
     * @brief Print locations of the qubit registers.
     */
    void print_reg_positions(std::ofstream& of) const;

    /**
     * @brief Save the qubit registers into the .circuit and .tex files.
     */
    void save_regs();

    /** @return pointer to the circuit. */
    YQCP get_the_circuit();

    /**
     * @brief Set the initial binary state of the circuit.
     */
    void set_init_binary_state(const bool& flag_mpi_bcast = false);

    void reset_init_vector(INIT_STATE__& state);

    /**
     * @brief Set initial amplitudes to specified qubits.
     * Elements \p ampl_vec_ in the state vector are filled starting from 
     * the elements corresponding to the low-priority qubit.
     * \p ampl_vec_real and \p ampl_vec_imag are assumed of the same size
     * @param[in] ampl_vec_real vector with real parts of amplitudes of size 2^nq;
     * @param[in] ampl_vec_imag vector with imaginary parts of amplitudes of size 2^nq;
     */
    void set_init_vector(YVQ ampl_vec_real, YVQ ampl_vec_imag);

    /**
     * @brief Set the qubit with id \p id_q to 1.
     * @param id_q id of the qubit to set to 1.
     */
    void set_qubit_state(YCU id_q);

    /**
     * @brief Set all qubits with id from \p ids_qs to 1.
     * @param ids_qs ids of the qubits to set to 1.
     */
    void set_qubit_state(YCVI ids_qs);

    /**
     * @brief Reset initial state of the circuit and remove all gates.
     */
    void reset();

    /**
     * @brief Reset Qureg. If the circuit is not yet allocated in memory, allocate it.
     */
    void reset_qureg();

    /**
     * @brief Empty saved binary initial states.
     */
    void empty_binary_states();

    /**
     * @brief Set zero state as an initial state of the circuit.
     */
    void prepare_zero_init_state();

    /**
     * @brief Set to 1 one of the qubits of a register. 
     * @param[in] name name of the register;
     * @param[in] id_reg_qubit qubit of the register to set to 1.
     */
    void set_reg_state(YCS name, YCI id_reg_qubit);

    /**
     * @brief Set to 1 several qubits of a register. 
     * @param[in] name name of the register;
     * @param ids_reg_qubits qubits of the register to set to 1.
     */
    void set_reg_state(YCS name, YCVI ids_reg_qubits);

    void read_structure_gate(
        YISS istr, YVI ids_target, qreal& par_gate, 
        YVI ids_unit, YVI ids_zero
    );
    void read_structure_gate(
        YISS istr, YVI ids_target, qreal& par_gate1, qreal& par_gate2,
        YVI ids_unit, YVI ids_zero
    );
    
    void read_structure_gate_adder_subtractor(YISS istr, YCS path_in, YCB flag_inv, YCI gate_type);
    void read_structure_gate_adder(YISS istr, YCS path_in, YCB flag_inv=false);
    void read_structure_gate_subtractor(YISS istr, YCS path_in, YCB flag_inv=false);
    void read_structure_gate_adder_qft(YISS istr, YCS path_in, YCB flag_inv=false);
    void read_structure_gate_subtractor_qft(YISS istr, YCS path_in, YCB flag_inv=false);
    void read_structure_gate_adder_fixed(YISS istr, YCS path_in, YCB flag_inv=false);
    void read_structure_gate_subtractor_fixed(YISS istr, YCS path_in, YCB flag_inv=false);
    void read_structure_gate_comparator_fixed(YISS istr, YCS path_in, YCB flag_inv=false);
    void read_structure_gate_swap(YISS istr, YCS path_in, YCB flag_inv=false);
    void read_structure_gate_fourier(YISS istr, YCS path_in, YCB flag_inv=false);
    void read_structure_sin(YISS istr, YCS path_in, YCB flag_inv=false);
    void read_structure_sinC(YISS istr, YCS path_in, YCB flag_inv=false);
    void read_structure_gate_phase_estimation(
        YISS istr, YCS path_in, std::map<std::string, YSQ>& ocs, YCB flag_inv
    );
    void read_structure_gate_qsvt(
        YISS istr, 
        std::map<std::string, YSQ>& ocs, 
        YCB flag_inv, 
        QuCF_complex_data& qucf_data,
        YCB flag_QETU = false
    );
    void read_structure_compression_gadget(
        YISS istr, 
        std::map<std::string, YSQ>& ocs, 
        YCB flag_inv,
        QuCF_complex_data& qucf_data
    );
    void read_structure_repeat(YISS istr, std::map<std::string, YSQ>& ocs, YCB flag_inv);
    void read_selector_power(YISS istr, std::map<std::string, YSQ>& ocs, YCB flag_inv=false);
    void read_structure_LCHS_QSP(YISS istr, std::map<std::string, YSQ>& ocs, YCB flag_inv=false);
    void read_structure_dirdec(YISS istr, YCB flag_inv=false);

    inline void read_global_control(YISS istr)
    {
        YVIv ids_unit, ids_zero;
        read_with_block(istr, ids_unit, ids_zero);
        blocks_ids_unit_.push_back(ids_unit);
        blocks_ids_zero_.push_back(ids_zero);
    }
    inline void remove_global_control()
    {
        blocks_ids_unit_.pop_back();
        blocks_ids_zero_.pop_back();
    }
    // inline void remove_repeat_block()
    // {
    //     Ns_repeat_.pop_back();
    //     remove_global_control();
    // }

    /**
     * @brief store indices of qubits found in a register pattern to the output vector \p ids_target.
     * ids_target[0] is the less significant qubit.
     * @param ids_target: qubits found in the register pattern;
     * @param flag_sort: if true, sort \p ids_target in non-descending order.
     */ 
    void read_reg_int(YISS istr, YVI ids_target, YCB flag_sort = true, YCS word_start=std::string());

    /**
     * @param ids_target: qubits found in the register pattern;
     * @param ids_target_e: all other qubits in the registers found in the pattern;
    */
    void read_reg_int(
        YISS istr, YVI ids_target, YVI ids_target_e,
        YCB flag_sort = true, YCS word_start=std::string()
    );

    /**
     * Read the rest of a gate description, where control nodes are described.
     * Add the control nodes to the corresponding vectors.
     * @param ids_unit: unit qubits; 
     * @param ids_zero: zero-control qubits;
     * @param id_element: 0 - gate, 1 - subcircuit, 2 - the with-structure;
    */
    void read_end_element(YISS istr, YVI ids_unit, YVI ids_zero, YCU id_element);
    inline void read_end_gate(YISS istr, YVI ids_unit, YVI ids_zero)
    {
        read_end_element(istr, ids_unit, ids_zero, 0);
    }
    inline void read_end_subcircuit(YISS istr, YVI ids_unit, YVI ids_zero)
    {
        read_end_element(istr, ids_unit, ids_zero, 1);
    }
    inline void read_with_block(YISS istr, YVI ids_unit, YVI ids_zero)
    {
        read_end_element(istr, ids_unit, ids_zero, 2);
    }
    // inline void read_repeat_block(YISS istr)
    // {
    //     std::string word;
    //     int N_repeat;
    //     istr >> word;
    //     N_repeat = get_value_from_word(word);
    //     Ns_repeat_.push_back(N_repeat);
    //     read_global_control(istr);
    // }

    template<class TGate>
    inline bool read_structure(YCS gate_name, YISS istr, YCB flag_inv=false)
    {
        YVIv ids_target, ids_unit, ids_zero;
        if(YMIX::compare_strings(gate_name, TGate::name_shared_, std::vector<std::string> {"X", "Y", "Z", "H"}))
        {
            qreal par_gate = nan("1");
            read_structure_gate(istr, ids_target, par_gate, ids_unit, ids_zero);
            for(auto const& id_target: ids_target) 
                add_sqg<TGate>(id_target, ids_unit, ids_zero, flag_inv);
            return true;
        }
        return false;
    }

    template<class TGate>
    inline bool read_structure(YCS gate_name, YISS istr, qreal& par_gate, YCB flag_inv=false)
    {
        YVIv ids_target, ids_unit, ids_zero;
        if(YMIX::compare_strings(gate_name, TGate::name_shared_, 
            std::vector<std::string> {"Rx", "Ry", "Rz", "Phase", "PhaseZero"}
        ))
        {
            read_structure_gate(istr, ids_target, par_gate, ids_unit, ids_zero);
            for(auto const& id_target: ids_target) 
                add_sq_rg<TGate>(id_target, par_gate, ids_unit, ids_zero, flag_inv);
            return true;
        }
        return false;
    }

    template<class TGate>
    inline bool read_structure(YCS gate_name, YISS istr, qreal& par_gate1, qreal& par_gate2, YCB flag_inv=false)
    {
        YVIv ids_target, ids_unit, ids_zero;
        if(YMIX::compare_strings(gate_name, TGate::name_shared_, std::vector<std::string> {"Rc"}))
        {
            read_structure_gate(istr, ids_target, par_gate1, par_gate2, ids_unit, ids_zero);
            for(auto const& id_target: ids_target) 
                add_sq_rg<TGate>(id_target, par_gate1, par_gate2, ids_unit, ids_zero, flag_inv);
            return true;
        }
        return false;
    }

    inline bool read_Rc_gadget(YCS gate_name, YISS istr, YCB flag_inv=false)
    {
        if(YMIX::compare_strings(gate_name, Rc_gadget__::name_shared_))
        {
            YVIv ids_target, ids_counter, ids_unit, ids_zero;
            std::string word;
            std::string filename;
            YVQv angles_ay, angles_az;

            // --- read target qubits ---
            read_reg_int(istr, ids_target);

            // --- read counter qubits ---
            read_reg_int(istr, ids_counter);

            // --- read the name of an .hdf5 file where angles of the gate are stored ---  
            istr >> word;
            filename = path_to_output_ + "/" + word + FORMAT_HDF5;

            // -- read the angles ---
            YMIX::H5File ff;
            ff.set_name(filename);
            ff.open_r();
            ff.read_vector(angles_ay, "ay", "angles");
            ff.read_vector(angles_az, "az", "angles");
            ff.close();

            // --- read the end of the gate structure description --- 
            read_end_gate(istr, ids_unit, ids_zero);

            // --- Create the gate ---
            for(auto const& id_target: ids_target)
                Rc_gadget(id_target, ids_counter, angles_ay, angles_az, ids_unit, ids_zero, flag_inv); 
            return true;
        }
        return false;
    }


    // add a single-qubit gate with several control gates:
    template<class TGate>
    YQCP add_sqg(YCI t, YCVI cs_unit = {}, YCVI cs_zero = {}, YCB flag_inv=false)
    {
        YSG oo = std::make_shared<TGate>(t);
        add_sq_core(oo, cs_unit, cs_zero, flag_inv);
        return get_the_circuit();
    }

    // add a single-qubit gate with one parameter with several control nodes: 
    template<class TGate>
    YQCP add_sq_rg(YCI t, YCQR a, YCVI cs_unit = {}, YCVI cs_zero = {}, YCB flag_inv=false)
    {
        YSG oo = std::make_shared<TGate>(t, a);
        add_sq_core(oo, cs_unit, cs_zero, flag_inv);
        return get_the_circuit();
    }

    // add a single-qubit gate with two parameters and with several control nodes: 
    template<class TGate>
    YQCP add_sq_rg(YCI t, YCQR a1, YCQR a2, YCVI cs_unit = {}, YCVI cs_zero = {}, YCB flag_inv=false)
    {
        YSG oo = std::make_shared<TGate>(t, a1, a2);
        add_sq_core(oo, cs_unit, cs_zero, flag_inv);
        return get_the_circuit();
    }


    /** Set Pauli X gate at \p t target qubit. 
     * @param[in] t target qubit;
     * @param[in] cs_unit 1-control qubits;
     * @param[in] cs_zero 0-control qubits;
     * @return pointer to the circuit.
     * */
    inline YQCP x(YCI t, YCVI cs_unit = {}, YCVI cs_zero = {}){ return add_sqg<X__>(t, cs_unit, cs_zero); }
    YQCP x(YCVI ts, YCVI cs_unit = {}, YCVI cs_zero = {});

    /** Set Pauli Y gate at \p t target qubit. */
    inline YQCP y(YCI t, YCVI cs_unit = {}, YCVI cs_zero = {}){ return add_sqg<Y__>(t, cs_unit, cs_zero); }
    YQCP y(YCVI ts, YCVI cs_unit = {}, YCVI cs_zero = {});

    /** Set Hadamard gate at \p t target qubit. */
    inline YQCP h(YCI t, YCVI cs_unit = {}, YCVI cs_zero = {}){ return add_sqg<H__>(t, cs_unit, cs_zero); }
    YQCP h(YCVI ts, YCVI cs_unit = {}, YCVI cs_zero = {});

    /** Set Pauli Z gate at \p t target qubit. */
    inline YQCP z(YCI t, YCVI cs_unit = {}, YCVI cs_zero = {}){ return add_sqg<Z__>(t, cs_unit, cs_zero); }
    YQCP z(YCVI ts, YCVI cs_unit = {}, YCVI cs_zero = {});

    /** Set a Rx-rotation gate.   
     * @param[in] t target qubit;
     * @param[in] a angle to rotate on;
     * @return pointer to the circuit.
     * */
    inline YQCP rx(YCI t, YCQR a, YCVI cs_unit = {}, YCVI cs_zero = {}, YCB flag_inv = false)
    { return add_sq_rg<Rx__>(t, a, cs_unit, cs_zero, flag_inv); }

    /** Set a Ry-rotation gate. */
    inline YQCP ry(YCI t, YCQR a, YCVI cs_unit = {}, YCVI cs_zero = {}, YCB flag_inv = false)
    { return add_sq_rg<Ry__>(t, a, cs_unit, cs_zero, flag_inv); }

    /** Set a Rz-rotation gate. */
    inline YQCP rz(YCI t, YCQR a, YCVI cs_unit = {}, YCVI cs_zero = {}, YCB flag_inv = false)
    { return add_sq_rg<Rz__>(t, a, cs_unit, cs_zero, flag_inv); }

    /** Set the Rc-rotation gate: Ry(ay).Rz(az)   
     * @param[in] t target qubit;
     * @param[in] az angle of the Rz-rotation;
     * @param[in] ay angle of the Ry-rotation;
     * @return pointer to the circuit.
     * */
    inline YQCP rc(YCI t, YCQR az, YCQR ay, YCVI cs_unit = {}, YCVI cs_zero = {}, YCB flag_inv = false)
    { return add_sq_rg<Rc__>(t, az, ay, cs_unit, cs_zero, flag_inv); }

    /** Set a phase shift gate.   
     * @param[in] t target qubit;
     * @param[in] a shift angle;
     * @return pointer to the circuit.
     * */
    inline YQCP phase(YCI t, YCQR a, YCVI cs_unit = {}, YCVI cs_zero = {}, YCB flag_inv = false)
    { return add_sq_rg<Phase__>(t, a, cs_unit, cs_zero, flag_inv); }

    inline YQCP phase_zero(YCI t, YCQR a, YCVI cs_unit = {}, YCVI cs_zero = {}, YCB flag_inv = false)
    { return add_sq_rg<PhaseZero__>(t, a, cs_unit, cs_zero, flag_inv); }


    inline YQCP Rc_gadget(
        YCI t, YCVI ids_counter, 
        YCVQ angles_ay, YCVQ angles_az, 
        YCVI cs_unit = {}, YCVI cs_zero = {}, 
        YCB flag_inv = false
    ){
        YSG oo = std::make_shared<Rc_gadget__>(t, angles_ay, angles_az, ids_counter);
        if(flag_inv) 
            oo->h_adjoint();
        oo->add_control_qubits(cs_unit, cs_zero);
        gates_.push_back(oo);
        // if(flag_layers_) 
        //     oo_layers_->add_gate(oo);
        return get_the_circuit();
    }


    /** @brief integer encoded to \p ts is incremented */
    YQCP adder_1(YCVI ts, YCVI cs_unit = {}, YCVI cs_zero = {}, YCB flag_inv= false);

    /** @brief integer encoded to \p ts + 2 */
    YQCP adder_2(YCVI ts, YCVI cs_unit = {}, YCVI cs_zero = {}, YCB flag_inv= false);

    /** @brief integer encoded to \p ts + 3 */
    YQCP adder_3(YCVI ts, YCVI cs_unit = {}, YCVI cs_zero = {}, YCB flag_inv= false);

    /** @brief integer encoded to \p ts is decremented */
    YQCP subtractor_1(YCVI ts, YCVI cs_unit = {}, YCVI cs_zero = {}, YCB flag_inv= false);

    /** @brief integer encoded to \p ts -2 */
    YQCP subtractor_2(YCVI ts, YCVI cs_unit = {}, YCVI cs_zero = {}, YCB flag_inv= false);

    /** @brief integer encoded to \p ts -3 */
    YQCP subtractor_3(YCVI ts, YCVI cs_unit = {}, YCVI cs_zero = {}, YCB flag_inv= false);

    /** @brief addition of two variables (v1 and v2) encoded to the registers \p ts1 and \p ts2;
     * the three registers, \p ts1, \p ts2 and \p ts3, must be of the same size;
     * the output sum (v1 + v2) is written to the qubits [ts2[:], ts3[0]].
     * The register \p ts3 must be initialized to the zero state.
     * The qubits \p ts3[1:] are used to store carry bits and are returned in the zero state.
     * @param flag_box if true, draw the operator as a box, not as a circuit;
     */
    YQCP adder(
        YCVI ts1, YCVI ts2, YCVI ts3, 
        YCVI cs_unit = {}, YCVI cs_zero = {},
        YCB flag_inv = false, YCB flag_box = false
    );


    /** @brief Subtraction of two variables (v1 and v2) encoded to the registers \p ts1 and \p ts2;
     * The three registers, \p ts1, \p ts2 and \p ts3, must be of the same size.
     * The output (v1 - v2) is written to the qubits [ts2[:], ts3[0]], where \p ts3[0] is the sign bit,
     * and |1>|00...0> corresponds to -1.
     * The register \p ts3 must be initialized to the zero state.
     * The qubits \p ts3[1:] are used to store carry bits and are returned in the zero state.
     */
    YQCP subtractor(
        YCVI ts1, YCVI ts2, YCVI ts3, 
        YCVI cs_unit = {}, YCVI cs_zero = {},
        YCB flag_inv = false
    );


    /** @brief addition of two variables (v1 and v2) encoded to the registers \p qs_v1 and \p qs_v2;
     * the two registers must be of the same size;
     * the output sum (v1 + v2) is written to the qubits [qs_v1[:], q_carry].
     * The qubit \p q_carry stores the carry bit of the sum.
     */
    YQCP adder_qft(
        YCVI qs_v1, YCVI qs_v2, YCI q_carry, 
        YCVI cs_unit = {}, YCVI cs_zero = {},
        YCB flag_inv = false, YCB flag_box = false
    );


    /** @brief subtraction of two variables (v1 and v2) encoded to the registers \p qs_v1 and \p qs_v2;
     * the two registers must be of the same size;
     * the output difference (v1 - v2) is written to the qubits [qs_v1[:], q_carry].
     * The qubit \p q_sign is inverted if the difference is negative.
     */
    YQCP subtractor_qft(
        YCVI qs_v1, YCVI qs_v2, YCI q_sign, 
        YCVI cs_unit = {}, YCVI cs_zero = {},
        YCB flag_inv = false, YCB flag_box = false
    );


    /** @brief addition of the unsigned integer \p int_sub to the unisgned integer encoded
     * into the qubits \p ids_target.
     * Results are stored back to \p ids_target.
     * The carry bit is written to the qubit \p id_carry.
     * The carry bit must be more significant than the qubits \p ids_target.
     */
    YQCP adder_fixed(
        YCVI ids_target, YCI id_carry, YCU int_sub, 
        YCVI cs_unit = {}, YCVI cs_zero = {},
        YCB flag_inv = false, YCB flag_box = false
    );


    /** @brief Integer, encoded within \p ids_target, - int_sub.
     * The results is written back to \p ids_target.
     * The sign bit is written to the qubit \p id_carry.
     * The carry bit must be more significant than the qubits \p ids_target.
     * Result |1>|00...0> means -1;
     * |1>|00...1>  -> -2 etc.
     */
    YQCP subtractor_fixed(
        YCVI ids_target, YCI id_carry, YCU int_sub, 
        YCVI cs_unit = {}, YCVI cs_zero = {},
        YCB flag_inv = false
    );


    /** @brief Compare unsigned integer \p int_sub with the unsigned integer encoded
     * within \p ids_target, call it uint_ref here.
     * If (\p int_sub > uint_ref), then the bit \p ids_carry[1] is inverted.
     * \p ids_carry must include two qubits, where \p ids_carry[1] is the most significant one.
     * The qubits \p ids_carry must be more significant than \p ids_target.
     * The states of the qubits \p ids_target and \p ids_carry[0] are remained untouched.
     */
    YQCP comparator_fixed(
        YCVI ids_target, YCVI ids_carry, YCU int_sub, 
        YCVI cs_unit = {}, YCVI cs_zero = {}, YCB flag_inv = false
    );


    /** @brief Add a swap operator between qubits \p t1 and \p t2. */
    YQCP swap(YCI t1, YCI t2, YCVI cs_unit = {}, YCVI cs_zero = {});

    /** @brief Quantum Fourier circuit placed to the qubits \p ts;
     * @param flag_box if true, draw the operator as a box, not as a circuit;
     */
    YQCP quantum_fourier(YCVI ts, YCVI cs_unit = {}, YCVI cs_zero = {}, YCB flag_inv = false, YCB flag_box = false);

    /** @brief Create sin(x), where x_j = alpha_0 + j*dx, j = [0, Nx-1], Nx = 2^size(\p conds),
     * dx = 2*alpha / Nx.
     * The sin(x) is created by Ry gates sitting on the ancilla \p anc[0] and controlled by qubits \p conds.
     * For the given integer j encoded as a bitstring in the qubits \p conds, 
     * the gate returns sin(x_j) as the amplitude of the state, 
     * where the ancilla \p anc[0] is in the zero state.
     */
    YQCP gate_sin(
        YCVI anc, 
        YCVI conds, 
        YCQR alpha_0, 
        YCQR alpha, 
        YCVI cs_unit = {}, YCVI cs_zero = {}, 
        YCB flag_inv = false, 
        YCB flag_box = false
    );

    /** @brief Create 
     *  sin(y) exp(i z), where 
     *      y_j = alpha_0_y + j*dy, dy = 2*alpha_y / N,  
     *      z_j = alpha_0_z + j*dz, dz = 2*alpha_z / N.          
     * Here, j = [0, N), N = 2^size(\p conds)
     */
    YQCP gate_sinC(
        YCVI anc, 
        YCVI conds, 
        YCQR alpha_0_y, YCQR alpha_y, 
        YCQR alpha_0_z, YCQR alpha_z,
        YCVI cs_unit = {}, YCVI cs_zero = {}, 
        YCB flag_inv = false
    );


    /** @brief Phase estimation (PE) operator to find an eigenphase of the operator \p A, 
     *         which sits on the qubits \p ta.
     *         The eigenstate of the circuit is prepared by the operator \p INIT.
     *         The final estimation is written to the qubits \p ty.
     * @param flag_box if true, draw the operator as a box, not as a circuit;
     */ 
    YQCP phase_estimation(
        YCVI ta, 
        YCCQ& A, 
        YCCQ& INIT,
        YCVI ty, 
        YCVI cs_unit = {}, YCVI cs_zero = {}, 
        YCB flag_inv = false,
        YCB flag_box = false
    );


    /** @brief Compression gadget.
     * @param[in] ids_counter qubits for the counter register;
     * @param[in] oc_U circuit which multiplication is to be computed;
     * @param[in] ids_U_target qubits where the operator \p oc_U sits;
     * @param[in] N_mult number of copies of \p oc_U in the product;
     * @param[in] flag_step_output whether state should be outputed after call to \p oc_U;
    */
    YQCP compression_gadget(
        GADGET_pars& data,
        YCVI ids_counter, 
        YCCQ& oc_U, 
        YCVI ids_U_target, 
        YCI N_mult, 
        YCB flag_step_output,
        YCVI cs_unit = {}, YCVI cs_zero = {},
        YCB flag_inv = false
    );

    /** @brief Repeat the oracle \p oc_U \p N_mult times.
     * @param[in] oc_U circuit which multiplication is to be computed;
     * @param[in] ids_U_target qubits where the operator \p oc_U sits;
     * @param[in] N_mult number of copies of \p oc_U in the product;
     * @param[in] flag_step_output whether state should be outputed after call to \p oc_U;
    */
    YQCP repeat(
        YCCQ& oc_U, 
        YCVI ids_U_target, 
        YCI N_mult, 
        YCVI cs_unit = {}, YCVI cs_zero = {},
        YCB flag_inv = false
    );

    /** @brief QSVT of the matrix encoded by the oracle \p BE, which sits on qubits \p qs_be.
     * The QSVT single rotations are placed at the qubit \p a_qsvt.
     * @param flag_box if true, draw the operator as a box, not as a circuit;
     */
    YQCP qsvt_def_parity(
        YCVQ phis,
        YCI a_qsvt,
        YCVI qs_be, 
        const std::shared_ptr<const QCircuit> BE,
        YCVI cs_unit = {}, YCVI cs_zero = {}, 
        YCB flag_inv = false,
        YCB flag_box = false
    );

    /** @brief QETU of the matrix encoded by the oracle \p BE, which sits on qubits \p qs_be.
     * The QSVT single rotations are placed at the qubit \p a_qsvt.
     * @param flag_box if true, draw the operator as a box, not as a circuit;
     */
    YQCP qetu_def_parity(
        YCVQ phis,
        YCI a_qsvt,
        YCVI qs_be, 
        const std::shared_ptr<const QCircuit> BE,
        YCVI cs_unit = {}, YCVI cs_zero = {}, 
        YCB flag_inv = false,
        YCB flag_box = false
    );


    /** @brief Hamiltonian simulation using the original (fully-coherent) QSP approach.
     * @param name_qsp_cricuit a unique name of the QSP circuit;
     * @param phis QSP angles for the whole exp(-i H t) function;
     * @param nt number of time intervals (n of repetitions of the QSP circuit)
     * @param a_qsp a qubit where the QSP rotations sit on;
     * @param a_qu a qubit for the qubitization;
     * @param qs_be qubits for the block-encoding oracle @p BE;
     * @param BE a block-encoding oracle
     * @param cs_unit qubits by which UNIT state the whole QSP circuit will be controlled.
     * @param cs_zero qubits by which ZERO state the whole QSP circuit will be controlled.
     * @param flag_inv if True, invert the circuit
    */
    YQCP qsp_ham(
        YCS name_qsp_cricuit,
        YCVQ phis,
        YCI nt,
        YCI a_qsp,
        YCI a_qu, 
        YCVI qs_be, 
        const std::shared_ptr<const QCircuit> BE,
        YCVI cs_unit = {}, YCVI cs_zero = {}, 
        YCB flag_inv = false
    );
    YQCP qsp_ham_opt2(
        YCS name_qsp_cricuit,
        YCVQ phis,
        YCI nt,
        YCI a_qsp,
        YCI a_qu, 
        YCVI qs_be, 
        const std::shared_ptr<const QCircuit> BE,
        YCVI cs_unit = {}, YCVI cs_zero = {}, 
        YCB flag_inv = false
    );

    /**
     * Selector |k><k> oc_U^k, where k = [0, 2^size(rs)).
     * @param rs selector qubits;
     * @param oc_U the unitary (subcircuit) whose powers will be computed;
     * @param ids_U_target qubits where the unitary sits;
    */
    YQCP selector_power(
        YCVI rs, 
        const std::shared_ptr<const QCircuit> oc_U, 
        YCVI ids_U_target,
        YCVI cs_unit, YCVI cs_zero,
        YCB flag_inv
    );

    /**
     * QSP-based Linear Combination of Hamiltonian simulations (LCHS_QSP).
     * @param Ut is the oracle simulating a short time interval.
     * @param ids_Ut are the target qubits of the oracle \p Ut.
     * @param Nt is the number of time intervals.
     * @param Ow is the oracle calculating the LCU weights.
     * @param ids_Ow are the target qubits of the oracle \p Ow.
    */
    YQCP LCHS_QSP(
        const std::shared_ptr<const QCircuit> Ut, YCVI ids_Ut, YCI Nt,
        const std::shared_ptr<const QCircuit> Ow, YCVI ids_Ow,
        YCVI cs_unit, YCVI cs_zero,
        YCB flag_inv
    );


    /**
     * Compute a profile by using multiple heavy controlled rotations Ry.
     * @param phis_y is the vector with rotation angles for the gate Ry;
     * @param ids_ctrl are the qubits controlling the rotations.
     * @param id_targ is the qubit storing the resulting profile.
    */
    YQCP DirDec_Y(
        YCVQ phis_y,
        YCVI ids_ctrl, 
        YCI id_targ,
        YCVI cs_unit, YCVI cs_zero,
        YCB flag_inv
    );



    /**
     * @brief Add a Stop gate to a quantum state at this point.
     * @param name is a name to identify this stop position.
     * @return YQCP 
     */
    inline YQCP add_stop_gate(YCS name)
    {
        YSG oo = std::make_shared<GStop__>(name);
        gates_.push_back(oo);
        return get_the_circuit();
    }

    inline void add_gate(std::shared_ptr<Gate__> oo){ gates_.push_back(oo); }

    /**
     * @brief Return the reference to the circuit statevector. 
     */
    void get_ref_to_state_vector(qreal*& state_real, qreal*& state_imag);

    /**
     * Return states with nonzero amplitudes.
     * @param[in] flag_ZeroHighPriorAnc if true, then assume that all ancillae are 
     * the high-priority qubits, and calculate 
     * only states where these ancillae are in the zero state.
     * @param[in] flag_zero_ampls if true, return zero amplitudes too.
     */
    void get_state(
        YMIX::StateVectorOut& out, 
        YCB flag_ZeroHighPriorAnc = false,
        YCB flag_zero_ampls = false
    );

    inline std::string get_name() const { return name_; }

    /**
     * @brief Get a number of ancilla qubits in the circuit. 
     */
    inline uint32_t get_na() const {return ancs_.size();}


    /**
     * Get names of nonancilla registers in the circuit.
     * @param reg_names_nonanc vector with names of nonancilla registers
     * ordered from the most to the least significant register
     * @param n_qubits is a vector with number of qubits in each nonancilla register:
     * order in the same way as @param reg_names_nonanc.
    */
    void get_nonancilla_regs(
        std::vector<std::string>& reg_names_nonanc,
        YVI n_qubits
    );

    inline int get_N_gates(){ return gates_.size(); }

private:
    qreal get_value_from_word(YCS word);
    void  qsvt_read_parameters(YCS filename, QSVT_pars& data);


    inline void add_sq_core(YSG& oo, YCVI cs_unit, YCVI cs_zero, YCB flag_inv)
    {
        if(flag_inv) 
            oo->h_adjoint();
        oo->add_control_qubits(cs_unit, cs_zero);
        gates_.push_back(oo);
        if(flag_layers_) 
            oo_layers_->add_gate(oo);
    }
    

    inline void get_id_qu_pattern(int& id_qu, YCS word, YCI nq_reg)
    {
        id_qu = get_value_from_word(word);
        if(id_qu < 0)
            id_qu = nq_reg + id_qu;

        if(id_qu >= nq_reg)
            throw std::string("\nqubit index, " + std::to_string(id_qu) + 
                ", is larger than the register size (" + std::to_string(nq_reg) + ")\n");
    }

    void read_reg_int_CORE(
        YISS istr, YVI ids_target, YVI ids_target_e,
        YCB flag_sort = true, YCS word_start=std::string(), YCB flag_e = false
    );

private:
    // if one addes a new property, do not forget to add it to the copy constructor.

    std::string name_; // name of the circuit;
    QuESTEnv env_; // QuEST execution environment;
    Qureg c_; // quantum register (circuit);
    YMIX::YTimer timer_; // timer for the circuit object;
    std::string path_to_output_; // path where output files to be written
    unsigned nq_; // number of qubits in the circuit;

    INIT_STATE__ init_state_; // initial state;

    std::string cfname_;  // name of the .circuit file;
    std::string texname_; // name of the .tex file;

    // registers:
    // regs_[rname][i] is the i-th qubit in the register "rname";
    // the 0-th qubit is the least significant in the register.
    std::map<std::string, YVIv> regs_; 

    // which registers are ancilla:
    std::map<std::string, bool> flags_anc_regs_;

    /** register names: the first name corresponds to the register at the top. */
    std::vector<std::string> regnames_; 

    // initial state of the circuit as an array of bits;
    // ib_state_[nq-1] corresponds to the least significant qubit 
    // (qubit at the bottom, the very right qubit)
    std::vector<short> ib_state_; 

    std::vector<YSG> gates_; // gates in the circuit;

    /** @brief Index of the starting gate to generate the circuit.*/
    int id_start_ = 0;

    std::vector<unsigned> standart_output_format_; // standart output format

    std::map<std::string, qreal> constants_;

    YVIv ancs_; // position of ancilla qubits (unordered vector);

    // object to organise gates in layers:
    std::shared_ptr<CircuitLayers__> oo_layers_;

    // circ_lines[i][j]: j-th phrase in the i-th line from the top:
    std::vector<std::vector<std::string>> tex_lines_;

    // names of qubts in the .tex:
    std::vector<std::string> tex_qubit_names_;

    // tex_noc_[i] = id of the next non-occupied position on the i-th qubit:
    std::vector<uint64_t> tex_noc_;

    // print or not .circuit files:
    bool flag_circuit_;

    //print or not .tex files:
    bool flag_tex_;

    // to calculate or not the layers;
    bool flag_layers_; 

    // if the circuit is allocated in memory:
    bool flag_circuit_allocated_;

    // // names of unique gates:
    // std::vector<std::string> unique_gates_names_;

    // blocks with global qubits for unit-control:
    std::vector<YVIv> blocks_ids_unit_;

    // blocks with global qubits for the zero-control:
    std::vector<YVIv> blocks_ids_zero_;

    // // Number of repeatitions for the "repeat" block:
    // std::vector<int> Ns_repeat_;

    bool flag_stop_gates_;
    bool flag_repeat_insert_;

};






#endif