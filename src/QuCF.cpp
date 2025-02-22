#include "../include/QuCF.h"
using namespace std;

QuCF__::QuCF__(
    const QuESTEnv& env, 
    YCS pname, 
    YCS path_to_inputs
) : BaseTool__(
    env, pname, path_to_inputs, 
    false,  // flag_compute_output (not used here),
    false,  // flag_print_output (not used here),
    false,  // default flag_circuit, 
    false,  // default flag_tex, 
    false,  // default flag_layers, 
    true    // default flag_hdf5
){
    format_file_ = FORMAT_ORACLE;
    sel_compute_output_   = "zero-ancillae";
    sel_print_output_     = "none";
    flag_init_state_file_ = false;
    flag_prob_ = false;
    flag_matrix_ = false;
    flag_zero_state_of_ = false;
    flag_output_qsp_ = false;
    flag_output_gadget_ = false;
    flag_stop_gates_ = true;
    flag_repeat_insert_ = false;
    flag_progress_bar_ = false;
    read_data();
}


QuCF__::~QuCF__()
{
    YMIX::print_log( "*** The oracletool has been destructed. ***");
}


void QuCF__::read_circuit_structure_from_file(YCS data)
{
    istringstream istr(data);
    string word;

    // read data
    while(istr >> word)
    {
        if(YMIX::compare_strings(word, "CONSTANTS"))
            read_constants(istr);

        if(YMIX::compare_strings(word, "OPTIONS"))
            read_options(istr);

        if(YMIX::compare_strings(word, "CIRCUITS_DECLARATION"))
            read_circuit_declaration(istr);

        if(YMIX::compare_strings(word, "CIRCUIT_STRUCTURE"))
            read_circuit_structure(istr);

        if(YMIX::compare_strings(word, "MAIN_CIRCUIT"))
            read_main_circuit(istr);
    }

    YMIX::print_log("Finish reading the .oracle file");
    YMIX::print_log("----------------------------------------------");

    // check if the circuit to launch is defined:
    if(!oc_to_launch_)
        throw "Error: The circuit to launch is not defined (the section MAIN_CIRCUIT is absent)."s;

    // write circuits to corresponding files:
    for(auto const& it: ocs_)
        it.second->print_gates();

    // set the output format for the circuit of interest:
    oc_to_launch_->set_standart_output_format();

    YMIX::print_log( "Oracle prepared.\n");
}


void QuCF__::read_constants(YISS istr)
{
    string word, constant_name;
    qreal constant_value;
    YMIX::print_log( "Reading constants...");
    try
    {
        while(istr >> word)
        {
            if(YMIX::compare_strings(word, "END_CONSTANTS"))
                break;

            // name of the constant
            constant_name = word; 

            // value of the constant
            istr >> word;
            constant_value = get_value_from_word(word);

            // save the constant:
            constants_[constant_name] = constant_value;
        }
    }
    catch(YCS e)
    {
        throw "Error in the CONSTANTS section when one reads the constant "s + constant_name + ": "s + e;
    }
}


void QuCF__::read_options(YISS istr)
{
    string word;
    YMIX::print_log( "Reading options...");

    try
    {
        while(istr >> word)
        {
            if(YMIX::compare_strings(word, "END_OPTIONS"))
                break;

            if(YMIX::compare_strings(word, "sel_compute_output"))
            {
                istr >> sel_compute_output_;
                if(!YMIX::compare_strings(sel_compute_output_, YVSv{"none", "all", "zero-ancillae"}))
                    throw string(
                        "unknown option ["s + sel_compute_output_ + "] for the selector sel_compute_output."s
                    );
                continue;
            }
            if(YMIX::compare_strings(word, "sel_print_output"))
            {
                istr >> sel_print_output_;
                if(!YMIX::compare_strings(sel_print_output_, YVSv{"none", "all", "zero-ancillae"}))
                    throw string(
                        "unknown option ["s + sel_print_output_ + "] for the selector sel_print_output."s
                    );

                if(
                    YMIX::compare_strings(sel_compute_output_, YVSv{"none", "all"}) &&
                    YMIX::compare_strings(sel_print_output_, YVSv{"zero-ancillae"})   
                )
                {
                    auto str_temp = YMIX::compare_strings(sel_compute_output_, "all") ? " separately": "";
                    cout << "\n-------------------------------------------------------------------------------\n"
                        << "WARNING: zero-ancillae states have not been computed" << str_temp << ":\n"  
                        << "sel_print_output is switched to " << sel_compute_output_ 
                        << "\n-------------------------------------------------------------------------------" << endl;
                    sel_print_output_ = sel_compute_output_;
                }

                if(
                    YMIX::compare_strings(sel_compute_output_, YVSv{"none", "zero-ancillae"}) &&
                    YMIX::compare_strings(sel_print_output_, YVSv{"all"})   
                )
                {
                    auto str_temp = 
                        YMIX::compare_strings(sel_compute_output_, "zero-ancillae") ?
                        "WARNING: not all states have been computed:\n":
                        "WARNING: states have not been computed:\n";
                    cout << "\n-------------------------------------------------------------------------------\n"
                        << str_temp  
                        << "sel_print_output is switched to " << sel_compute_output_ 
                        << "\n-------------------------------------------------------------------------------" << endl;
                    sel_print_output_ = sel_compute_output_;
                }

                continue;
            }

            if(YMIX::compare_strings(word, "flag_output_qsp"))
            {
                istr >> flag_output_qsp_;
                continue;
            }

            if(YMIX::compare_strings(word, "flag_output_gadget"))
            {
                istr >> flag_output_gadget_;
                continue;
            }

            if(YMIX::compare_strings(word, "flag_circuit"))
            {
                istr >> flag_circuit_;
                continue;
            }

            if(YMIX::compare_strings(word, "flag_tex"))
            {
                istr >> flag_tex_;
                continue;
            }

            if(YMIX::compare_strings(word, "tex_CL"))
            {
                istr >> YGV::tex_circuit_length;
                continue;
            }

            if(YMIX::compare_strings(word, "flag_matrix"))
            {
                istr >> flag_matrix_;
                continue;
            }

            if(YMIX::compare_strings(word, "flag_stop_gates"))
            {
                istr >> flag_stop_gates_;
                continue;
            }

            if(YMIX::compare_strings(word, "flag_repeat_insert"))
            {
                istr >> flag_repeat_insert_;
                continue;
            }

            if(YMIX::compare_strings(word, "flag_progress_bar"))
            {
                istr >> flag_progress_bar_;
                continue;
            }
        }

        // correct the printing options:
        string sel_comp = "none";
        if(YMIX::compare_strings(sel_compute_output_, sel_comp))
            if(!YMIX::compare_strings(sel_print_output_, sel_comp))
                sel_print_output_ = sel_comp;

        sel_comp = "zero-ancillae";
        if(YMIX::compare_strings(sel_compute_output_, sel_comp))
            if(YMIX::compare_strings(sel_print_output_, "all"))
                sel_print_output_ = sel_comp; // "all" -> "zero-ancillae"

        // output some options:
        YMIX::print_log( "\n--- Initial OPTIONS ---");
        if(YMIX::compare_strings(sel_compute_output_, "none"))
            YMIX::print_log("-> do not compute output states.");
        if(YMIX::compare_strings(sel_compute_output_, "all"))
            YMIX::print_log("-> compute all output states.");
        if(YMIX::compare_strings(sel_compute_output_, "zero-ancillae"))
            YMIX::print_log("-> compute only output states, where all ancillae are in the zero state.");
        if(!flag_circuit_)        YMIX::print_log("-> do not write the " + FORMAT_CIRCUIT + " files.");
        if(!flag_tex_)            YMIX::print_log("-> do not write the " + FORMAT_TEX + " file.");
        if(flag_tex_) 
            YMIX::print_log(
                "-> length of a single row in the .tex file = " + to_string(YGV::tex_circuit_length)
            );
        if(flag_stop_gates_)
            YMIX::print_log("-> Include stop gates.");
        else
            YMIX::print_log("-> Do not include stop gates.");
        if(flag_repeat_insert_)
            YMIX::print_log("-> In the repeat-gate, just INSERT gates.");
        YMIX::print_log( "---\n");
    }
    catch(YCS e)
    {
        throw "Error in the OPTIONS section: "s + e;
    }
}


qreal QuCF__::get_value_from_word(YCS word)
{
    if(word.find("<") == string::npos)
    {
        istringstream sstr(word);
        qreal res_value;
        if(!(sstr >> res_value))
            throw "wrong format"s;
        return res_value;
    }

    unsigned first = word.find("<");
    unsigned last = word.find(">");
    string const_name = word.substr(first+1,last-first-1);

    if(constants_.find(const_name) == constants_.end())
        throw "The constant with the name ["s + const_name + "] is not found."s;

    return constants_[const_name];
}


void QuCF__::read_circuit_declaration(YISS istr)
{
    string word, circ_name, current_field;
    int n_regs, n_qubits_in_reg;

    YMIX::print_log("Reading declaration of circuits...", 0, true);

    try
    {
        while(istr >> word)
        {
            if(YMIX::compare_strings(word, "END_CIRCUITS_DECLARATION"))
                break;

            // name of the circuit
            circ_name = word; 
            
            // number of registers in the circuit
            current_field = "the number of registers"s; 
            istr >> word;
            n_regs = get_value_from_word(word);

            // registers are organized from the high to low priority
            // (from top to bottom qubits):
            vector<int> saved_nqs_in_regs(n_regs);
            vector<string> saved_reg_names(n_regs);
            vector<bool> anc_flags(n_regs);
            int nq_circ = 0;
            YVIv ancs;
            for(unsigned i = 0; i < n_regs; i++)
            {
                string str_value = to_string(i) + "-th register";
                if(i == 1)
                    str_value = to_string(i) + "-st register";
                if(i == 2)
                    str_value = to_string(i) + "-nd register";

                // register name
                current_field = "the name of the " + str_value;
                istr >> saved_reg_names[i];

                // number of qubits in the register
                current_field = "the number of qubits in the " + str_value;
                istr >> word;
                saved_nqs_in_regs[i] = get_value_from_word(word);

                nq_circ += saved_nqs_in_regs[i];

                // ancilla flag
                current_field = "the ancilla flag of the " + str_value;
                istr >> word;
                anc_flags[i] = int(get_value_from_word(word));
            }

            // create a circuit if necessary; choose a circuit:
            if(ocs_.find(circ_name) != ocs_.end())
            {
                ostringstream inf;
                inf << "Warning: The circuit " << circ_name <<  " is declared several times. " <<
                        "The first declaration is taken.\n";
                YMIX::print_log( inf.str());          
            }
            else
            {
                // cout << "HERE 2: circuit [" << circ_name << "]: " << flag_tex_ << endl;

                // at this moment, circuits are NOT allocated in memory:
                ocs_[circ_name] = make_shared<QCircuit>(
                    circ_name, env_, path_inputs_, nq_circ, 
                    constants_,
                    flag_circuit_, flag_tex_, flag_layers_, 
                    flag_stop_gates_, flag_repeat_insert_, flag_progress_bar_
                );
            }
                
            YSQ oc = ocs_[circ_name];

            // add the registers to the chosen circuit
            for(unsigned i = 0; i < n_regs; i++)
                oc->add_register(saved_reg_names[i], saved_nqs_in_regs[i], anc_flags[i]);
            oc->save_regs();
        }
    }
    catch(YCS e)
    {
        throw "Circuit Declaration: Error when one reads "s + current_field + " in the circuit "s + circ_name + ": "s + e;
    }

    // information about the declared circuits
    for(const auto& it: ocs_) 
    {
        ostringstream inf;
        inf << it.first << " with " << it.second->get_n_qubits() << " qubit(s):";
        YMIX::print_log(inf.str(), 1, true); 

        inf.str(""); inf.clear();
        auto reg_names = it.second->get_reg_names();
        inf << "with " << reg_names.size() << " register(s):\n";
        for(auto const& reg_name: reg_names)
            inf << "register " << reg_name << " with " << 
                it.second->get_nq_in_reg(reg_name) << " qubit(s);\n";
        YMIX::print_log(inf.str(), 2, true); 
    }
}


void QuCF__::read_circuit_structure(YISS istr, YSQ* oc_ext)
{
    string word;
    string curr_circuit_name;
    int n_regs, n_qubits_in_reg;
    YMIX::YTimer timer;

    // choose the circuit according to its name:
    YSQ oc;
    if(oc_ext == nullptr)
    {
        istr >> word;
        YMIX::print_log("----> Reading structure of the circuit " + word, 0, false);
        if(ocs_.find(word) == ocs_.end()) 
        {
            YMIX::print_log(
                "Warning: No circuit with the name " + word + " has been declared. Skip it.", 1
            );
            return;
        }
        curr_circuit_name = word;
        // YSQ oc = ocs_[curr_circuit_name];
        oc = ocs_[curr_circuit_name];
    }
    else
    {
        oc = *oc_ext;
    }  
    timer.StartPrint("reading...", 0, true, false);

    // --- gates in the circuit ---
    while(istr >> word)
    {
        // cout << "word: " << word << endl;

        if(YMIX::compare_strings(word, "END_CIRCUIT_STRUCTURE"))
            break;

        if(YMIX::compare_strings(word, "gate"))
            read_gate(istr, oc);

        if(YMIX::compare_strings(word, "igate"))
            read_gate(istr, oc, true);

        if(YMIX::compare_strings(word, "circuit"))
            read_subcircuit(istr, oc, false);

        if(YMIX::compare_strings(word, "icircuit"))
            read_subcircuit(istr, oc, true);

        if(YMIX::compare_strings(word, "file"))
            read_gates_from_file(istr, oc);

        if(YMIX::compare_strings(word, "with"))
            oc->read_global_control(istr);

        if(YMIX::compare_strings(word, "end_with"))
            oc->remove_global_control();

        // if(YMIX::compare_strings(word, "repeat"))
        //     oc->read_repeat_block(istr);

        // if(YMIX::compare_strings(word, "repeat"))
        //     oc->remove_repeat_block();
    }
    timer.StopPrint();
    YMIX::print_log("--------------------------------------------", 0, false);
}


void QuCF__::read_gate(YISS istr, YPQC oc, YCB flag_inv)
{
    string gate_name;
    qreal par_gate, par_gate2;
     
    istr >> gate_name;

    // cout << "here: " << gate_name << endl;

    try
    {
        if(oc->read_structure<X__>(gate_name, istr, flag_inv)) return;
        if(oc->read_structure<Y__>(gate_name, istr, flag_inv)) return;
        if(oc->read_structure<Z__>(gate_name, istr, flag_inv)) return;
        if(oc->read_structure<H__>(gate_name, istr, flag_inv)) return;

        if(oc->read_structure<Rx__>(gate_name, istr, par_gate, flag_inv)) return;
        if(oc->read_structure<Ry__>(gate_name, istr, par_gate, flag_inv)) return;
        if(oc->read_structure<Rz__>(gate_name, istr, par_gate, flag_inv)) return;
        if(oc->read_structure<Phase__>(gate_name, istr, par_gate, flag_inv)) return;

        if(oc->read_structure<Rc__>(gate_name, istr, par_gate, par_gate2, flag_inv)) return;

        if(oc->read_Rc_gadget(gate_name, istr, flag_inv)) return;

        if(YMIX::compare_strings(gate_name, YVSv{"incrementor", "adder1"}))
        {
            oc->read_structure_gate_adder_subtractor(istr, path_inputs_, flag_inv, 1);
            return;
        }
        if(YMIX::compare_strings(gate_name, "adder2"))
        {
            oc->read_structure_gate_adder_subtractor(istr, path_inputs_, flag_inv, 2);
            return;
        }
        if(YMIX::compare_strings(gate_name, "adder3"))
        {
            oc->read_structure_gate_adder_subtractor(istr, path_inputs_, flag_inv, 3);
            return;
        }
        if(YMIX::compare_strings(gate_name, YVSv{"decrementor", "subtractor1"}))
        {
            oc->read_structure_gate_adder_subtractor(istr, path_inputs_, flag_inv, -1);
            return;
        }
        if(YMIX::compare_strings(gate_name, "subtractor2"))
        {
            oc->read_structure_gate_adder_subtractor(istr, path_inputs_, flag_inv, -2);
            return;
        }
        if(YMIX::compare_strings(gate_name, "subtractor3"))
        {
            oc->read_structure_gate_adder_subtractor(istr, path_inputs_, flag_inv, -3);
            return;
        }
        if(YMIX::compare_strings(gate_name, "adder"))
        {
            oc->read_structure_gate_adder(istr, path_inputs_, flag_inv);
            return;
        }
        if(YMIX::compare_strings(gate_name, "subtractor"))
        {
            oc->read_structure_gate_subtractor(istr, path_inputs_, flag_inv);
            return;
        }
        if(YMIX::compare_strings(gate_name, "adder_qft"))
        {
            oc->read_structure_gate_adder_qft(istr, path_inputs_, flag_inv);
            return;
        }
        if(YMIX::compare_strings(gate_name, "subtractor_qft"))
        {
            oc->read_structure_gate_subtractor_qft(istr, path_inputs_, flag_inv);
            return;
        }
        if(YMIX::compare_strings(gate_name, "AdderFixed"))
        {
            oc->read_structure_gate_adder_fixed(istr, path_inputs_, flag_inv);
            return;
        }
        if(YMIX::compare_strings(gate_name, "SubtractorFixed"))
        {
            oc->read_structure_gate_subtractor_fixed(istr, path_inputs_, flag_inv);
            return;
        }
        if(YMIX::compare_strings(gate_name, "ComparatorFixed"))
        {
            oc->read_structure_gate_comparator_fixed(istr, path_inputs_, flag_inv);
            return;
        }
        if(YMIX::compare_strings(gate_name, "swap"))
        {
            oc->read_structure_gate_swap(istr, path_inputs_, flag_inv);
            return;
        }
        if(YMIX::compare_strings(gate_name, "Fourier"))
        {
            oc->read_structure_gate_fourier(istr, path_inputs_, flag_inv);
        }
        if(YMIX::compare_strings(gate_name, "sin"))
        {
            oc->read_structure_sin(istr, path_inputs_, flag_inv);
        }
        if(YMIX::compare_strings(gate_name, "sinCO"))
        {
            oc->read_structure_sinC(istr, path_inputs_, flag_inv);
        }
        // if(YMIX::compare_strings(gate_name, "sinCTRL"))
        // {
        //     oc->read_structure_sinCTRL(istr, path_inputs_, flag_inv);
        // }
        if(YMIX::compare_strings(gate_name, "PE"))
        {
            oc->read_structure_gate_phase_estimation(istr, path_inputs_, ocs_, flag_inv);
        }
        // if(YMIX::compare_strings(gate_name, "PE_BE"))
        // {
        //     oc->read_structure_gate_phase_estimation_BE(istr, path_inputs_, ocs_, flag_inv);
        // }
        if(YMIX::compare_strings(gate_name, "QSVT"))
        {
            oc->read_structure_gate_qsvt(istr, ocs_, flag_inv, cd_);
        }
        if(YMIX::compare_strings(gate_name, "QETU"))
        {
            bool flag_QETU = true;
            oc->read_structure_gate_qsvt(istr, ocs_, flag_inv, cd_, flag_QETU);
        }
        if(YMIX::compare_strings(gate_name, "CompressionGadget"))
        {
            oc->read_structure_compression_gadget(istr, ocs_, flag_inv, cd_);
        }
        if(YMIX::compare_strings(gate_name, "repeat"))
        {
            oc->read_structure_repeat(istr, ocs_, flag_inv);
        }
        if(YMIX::compare_strings(gate_name, "SelectorPower"))
        {
            oc->read_selector_power(istr, ocs_, flag_inv);
        }
        if(YMIX::compare_strings(gate_name, "LCHS_QSP"))
        {
            oc->read_structure_LCHS_QSP(istr, ocs_, flag_inv);
        }
        if(YMIX::compare_strings(gate_name, "DirDec"))
        {
            oc->read_structure_dirdec(istr, flag_inv);
        }
    }
    catch(YCS e)
    {
        ostringstream ostr;
        ostr << "--- Error in the structure of the circuit " << oc->get_name() 
             << ", in the gate " << gate_name << " ---\n" <<
                "" << e << " ---\n";
        throw ostr.str();
    }
}


void QuCF__::read_subcircuit(YISS istr, YPQC oc, YCB flag_inv)
{
    string subcircuit_name, word;
    string curr_circuit_name = oc->get_name();
    int id_qubit;
    unsigned nq_sub;
    bool flag_skip = false;
    bool flag_plain = false;

    istr >> subcircuit_name;
    if(ocs_.find(subcircuit_name) == ocs_.end())
    {
        string warn_line;
        warn_line = "\n\n-------------------------------------------------------------------------------\n";
        warn_line += "--- Warning: setting the structure of the circuit [" + curr_circuit_name + "] ---\n";
        warn_line += "The subcircuit [" + subcircuit_name + "] is not found. We skip it.\n";
        warn_line += "-------------------------------------------------------------------------------\n";
        YMIX::print_log( warn_line);
        flag_skip = true;
    }
    if(YMIX::compare_strings(subcircuit_name, curr_circuit_name))
    {
        string warn_line;
        warn_line = "\n\n-------------------------------------------------------------------------------\n";
        warn_line += "--- Warning: setting the structure of the circuit [" + curr_circuit_name + "] ---\n";
        warn_line += "Self-insertion: circuit cannot include itself as a subcircuit. We skip it.\n";
        warn_line += "-------------------------------------------------------------------------------\n";
        YMIX::print_log( warn_line);
        flag_skip = true;
    }
    if(ocs_[subcircuit_name]->get_N_gates() == 0)
    {
        string warn_line;
        warn_line = "\n\n-------------------------------------------------------------------------------\n";
        warn_line += "--- Warning: setting the structure of the circuit [" + curr_circuit_name + "] ---\n";
        warn_line += "The subcircuit [" + subcircuit_name + "] is empty. We skip it.\n";
        warn_line += "-------------------------------------------------------------------------------\n";
        YMIX::print_log( warn_line);
        flag_skip = true;
    }

    if(!flag_skip)
    {
        YSQ oc_sub  = ocs_[subcircuit_name];
        int plain_int;

        istr >> word;

        // check whether subcircuit qubits are connected to 
        //         the current circuit qubits in a trivial way:
        try 
        {
            plain_int = stoi(word);
            if(plain_int == -1) flag_plain = true;
        }
        catch(...)
        {   
            flag_plain = false;
        }

        nq_sub = oc_sub->get_n_qubits();
        vector<int> ids_q(nq_sub);
        if(flag_plain)
            for(unsigned i = 0; i < nq_sub; i++)
                ids_q[i] = i;
        else
        {
            ids_q = YVIv {};
            oc->read_reg_int(istr, ids_q, true, word);
            if(ids_q.size() != nq_sub)
            {
                string err_line;
                err_line  = "--- Error: setting the structure of the circuit " + curr_circuit_name + " ---\n";
                err_line += "The subcircuit " + subcircuit_name + " has " + to_string(nq_sub) + " qubits, while" + 
                    " one has indicated " + to_string(ids_q.size()) + " qubits to connect to.";
                throw err_line;
            }
        }

        // read the end of the subcircuit (control qubits) and find the end_circuit keyword:
        YVIv ids_unit, ids_zero;
        oc->read_end_subcircuit(istr, ids_unit, ids_zero);

        // add gates from the subcircuit to the parent circuit:
        oc->copy_gates_from(oc_sub, ids_q, YSB(nullptr), ids_unit, ids_zero, flag_inv);
    }
}


void QuCF__::read_gates_from_file(YISS istr, YPQC oc)
{
    string data;
    string file_name;
    istr >> file_name;
    file_name += ".oracle";
    YMIX::read_input_file(data, file_name);

    istringstream istr_file(data);
    read_circuit_structure(istr_file, &oc);
}


void QuCF__::read_main_circuit(YISS istr)
{
    YMIX::print_log("Reading MAIN_CIRCUIT section...", 0, true);
    string word, name_of_main_circuit;
    bool flag_init_file = false;

    // define a circuit to launch:
    istr >> name_of_main_circuit;
    if(ocs_.find(name_of_main_circuit) == ocs_.end())
        throw "Error: The circuit with a name \"" + name_of_main_circuit + "\" is not found.\n" + 
            "This circuit cannot be set as a circuit to launch.";
    oc_to_launch_ = ocs_[name_of_main_circuit];

    try
    {
        // to read the initial state directly from the .oracle file
        while(istr >> word)
        {
            // to read the initial state from the .init file:
            if(YMIX::compare_strings(word, "INIT_FILE"))
            {
                read_state_init_file();
                flag_init_file = true;
                continue;
            }
                
            if(YMIX::compare_strings(word, "INPUT_STATE"))
            {
                if(!flag_init_file)
                    read_state(istr);
            }
       
            if(YMIX::compare_strings(word, "compute_prob"))
            {
                istr >> flag_prob_;
                oc_to_launch_->read_reg_int(istr, focus_qubits_);

                stringstream sstr;
                copy(
                    focus_qubits_.begin(), 
                    focus_qubits_.end(), 
                    std::ostream_iterator<int>(sstr, " ")
                );
                YMIX::print_log("-> compute probabilites of the states on the following qubits: " + sstr.str());
                continue;
            }

            if(YMIX::compare_strings(word, "output_zero_state_of"))
            {
                flag_zero_state_of_ = true;
                YVIv ids_qs;
                oc_to_launch_->read_reg_int(istr, ids_qs);

                int nq = oc_to_launch_->get_n_qubits(); 
                qs_chosen_state_.resize(nq);
                std::fill(qs_chosen_state_.begin(), qs_chosen_state_.end(), -1);
                for(auto id_qu: ids_qs)
                    qs_chosen_state_[nq - id_qu - 1] = 0;

                stringstream sstr;
                copy(
                    qs_chosen_state_.begin(), 
                    qs_chosen_state_.end(), 
                    std::ostream_iterator<int>(sstr, " ")
                );
                YMIX::print_log("-> compute probabilites of the states on the following qubits: " + sstr.str());
                continue;
            }

            if(YMIX::compare_strings(word, "END_MAIN_CIRCUIT"))
                break;
        }

        if(!flag_init_file && init_states_.size() == 0)
            YMIX::print_log("\n<<<WARNING: there are not input states.>>>\n");

    }
    catch(YCS e)
    {
        ostringstream ostr;
        ostr << "--- Error while reading initial states: ---"  <<
                "\n--- " << e << " ---\n";
        throw ostr.str();
    }
}


void QuCF__::read_state(YISS istr)
{
    string word;
    std::map<string, vector<int>> one_state;
    vector<int> ids_qs;

    oc_to_launch_->read_reg_int(istr, ids_qs);
    one_state[YGV::reg_whole_circuit] = ids_qs;
    init_states_.push_back(one_state); 
}


void QuCF__::read_state_init_file()
{
    int N;
    flag_init_state_file_ = true;

    string fname_init = path_inputs_ + "/" + pname_ + FORMAT_INIT;
    cout << "Reading the file: " << path_inputs_ + "/" + pname_ + FORMAT_INIT << endl;
    YMIX::read_init_state(fname_init, init_ampl_vec_real_, init_ampl_vec_imag_);
}


void QuCF__::launch()
{
    YMIX::print_log( 
        "\n-------------------------------------------------------------------------------"s + 
        "\n-------------------------------------------------------------------------------"s
    );
    YMIX::print_log( 
        "--- Analysis of the circuit [" + oc_to_launch_->get_name() + "] ---"
    );
    
    // working circuit object:
    shared_ptr<QCircuit> u_work = oc_to_launch_;
    
    // ---- store basic data:
    hfo_.open_w();

    // number of qubits
    hfo_.add_scalar(u_work->get_n_qubits(), "nq", "basic");
    hfo_.add_scalar(u_work->get_na(), "na", "basic");
    YMIX::print_log("nq: " + to_string(u_work->get_n_qubits()));
    YMIX::print_log("na: " + to_string(u_work->get_na()));

    // number of gates
    hfo_.add_scalar(u_work->get_N_gates(), "N_gates", "basic");
    YMIX::print_log("N-gates: " + to_string(u_work->get_N_gates()));

    // register names
    string res_lin = "";
    for(auto const& reg_name: u_work->get_reg_names())
        res_lin += reg_name + ", ";
    res_lin.pop_back(); res_lin.pop_back();
    hfo_.add_scalar(res_lin, "register-names", "basic");

    // number of qubits in every register:
    hfo_.add_vector(u_work->get_standart_output_format(), "register-nq", "basic");

    // if QSVT, store its parameters:
    auto& qsvt_d = cd_.qsvt;
    if(!qsvt_d.empty())
    {
        YMIX::print_log("Saving QSVT parameters...");
        hfo_.add_group("qsvt");
        hfo_.add_scalar(qsvt_d.size(), "n-qsvt-circuits", "qsvt");

        int counter_qsvt = -1;
        for(auto const& [name_qsvt_gate, qsvt_data_one] : qsvt_d)
        {
            ++counter_qsvt;
            string name_gr = name_qsvt_gate;

            hfo_.add_scalar(name_gr, "name-"s + to_string(counter_qsvt), "qsvt");

            hfo_.add_group(name_gr);
            hfo_.add_scalar(name_qsvt_gate,         "name",   name_gr);
            hfo_.add_scalar(qsvt_data_one.type,     "type",   name_gr);
            hfo_.add_scalar(qsvt_data_one.eps_qsvt, "eps",    name_gr);
            hfo_.add_scalar(qsvt_data_one.parity,   "parity", name_gr);
            hfo_.add_scalar(qsvt_data_one.rescaling_factor,   "rescaling_factor", name_gr);
            hfo_.add_scalar(qsvt_data_one.f_par, "function-parameter", name_gr);

            // if(qsvt_data_one.parity == 0)
            // {
            //     hfo_.add_scalar(qsvt_data_one.angles_phis_even, "angles-even", name_gr);
            // }
            // if(qsvt_data_one.parity == 1)
            // {
            //     hfo_.add_scalar(qsvt_data_one.angles_phis_odd, "angles-odd", name_gr);
            // }

            // if(YMIX::compare_strings(qsvt_data_one.type, "matrix-inversion"))
            // {
            //     hfo_.add_scalar(qsvt_data_one.f_par, "kappa", name_gr);
            //     hfo_.add_scalar(qsvt_data_one.angles_phis_odd, "angles-odd", name_gr);
            // }
            // if(YMIX::compare_strings(qsvt_data_one.type, "xgaussian"))
            // {
            //     hfo_.add_scalar(qsvt_data_one.f_par, "mu", name_gr);
            //     hfo_.add_scalar(qsvt_data_one.angles_phis_odd, "angles-odd", name_gr);
            // }
            // if(YMIX::compare_strings(qsvt_data_one.type, "gaussian-arcsin"))
            // {
            //     hfo_.add_scalar(qsvt_data_one.f_par, "mu", name_gr);
            //     hfo_.add_scalar(qsvt_data_one.angles_phis_even, "angles-even", name_gr);
            // }
            // if(YMIX::compare_strings(qsvt_data_one.type, "QSVT-ham"))
            // {
            //     hfo_.add_scalar(qsvt_data_one.f_par, "dt", name_gr);
            //     hfo_.add_scalar(qsvt_data_one.n_repeat, "nt", name_gr);
            //     hfo_.add_scalar(qsvt_data_one.angles_phis_odd, "angles-odd", name_gr);
            //     hfo_.add_scalar(qsvt_data_one.angles_phis_even, "angles-even", name_gr);
            // }
            if(YMIX::compare_strings(qsvt_data_one.type, "QSP-ham"))
            {
                hfo_.add_scalar(qsvt_data_one.f_par, "dt", name_gr);
                hfo_.add_scalar(qsvt_data_one.n_repeat, "nt", name_gr);
                // hfo_.add_scalar(qsvt_data_one.angles_phis_arbitrary, "angles", name_gr);
            }
        }
    }

    // Save data of gadgets:
    auto& gadgets_d = cd_.gadgets;
    if(!gadgets_d.empty())
    {
        YMIX::print_log("Saving Gadgets' parameters...");
        for(auto const& [gadget_name, gadget_data] : gadgets_d)
        {
            hfo_.add_group(gadget_name);
            hfo_.add_scalar(gadget_data.type, "type", gadget_name);
            if(YMIX::compare_strings(gadget_data.type, "compression"))
            {
                hfo_.add_scalar(gadget_data.N_mult, "N_mult", gadget_name);
                hfo_.add_vector(gadget_data.counter_qubits, "counter-qubits", gadget_name);
            }
        }
    }

    // number of initial states:
    uint32_t n_init_states = flag_init_state_file_ ? 1: init_states_.size();
    hfo_.add_scalar(n_init_states, "n-init-states", "states");

    // constants:
    for(auto const& [key, val] : constants_)
    {
        hfo_.add_scalar(val, key, "constants");
    }

    hfo_.close();

    // -------------------------------------------------------------------------------------------
    // --- Matrix construction ---
    if(flag_matrix_)
        calc_matrix(u_work, false);
    
    // -------------------------------------------------------------------------------------------
    // --- Analyse and Store if necessary initial and output states of the circuit ---
    int count_init_state = 0;
    if(flag_init_state_file_)
    {
        // remove gates from the previous launch:
        u_work->reset_qureg();

        cout << "size of the initial vector is " << init_ampl_vec_real_.size() << endl;
        u_work->set_init_vector(init_ampl_vec_real_, init_ampl_vec_imag_);
        calc(u_work, count_init_state);
    }
    else
    {
        for(auto const& state: init_states_)
        {
            // empty previous initial binary states:
            u_work->empty_binary_states();

            // remove gates from the previous launch:
            u_work->reset_qureg();
            
            // set initial states:
            for(auto const& reg: state)
                u_work->set_reg_state(reg.first, reg.second);
            u_work->set_init_binary_state();

            // print input and output states:
            calc(u_work, count_init_state);
            ++count_init_state;
        }
    }
}


void QuCF__::calc(shared_ptr<QCircuit>& u_work, YCI count_init_state)
{
    YMIX::YTimer timer_comp;
    // --- Print the initial state ---
    {
        YMIX::StateVectorOut outF;
        u_work->get_state(outF);
        YMIX::print_log(
            "..........................................................................\n" + 
            "...Initial state "s + to_string(count_init_state) + "...\n" + outF.str_wv
        );
    
        // --- Store the initial state ---
        hfo_.open_w();
        hfo_.add_vector(outF.ampls,  "initial-amplitudes-"s + to_string(count_init_state), "states");
        hfo_.add_matrix(outF.states, "initial-states-"s + to_string(count_init_state),     "states");
        hfo_.close(); 
    }

    // --- Print output states ---
    if(!YMIX::compare_strings(sel_compute_output_, "none") or flag_output_qsp_ or flag_output_gadget_)
    {
        uint64_t id_current_gate = 0;
        string stop_point_name;
        YMIX::StateVectorOut outF, outZ, outZ_chosen, outZ_qsp, outZ_gadget;

        map<string, int> counters_t;
        for(auto const& [name_qsvt_gate, qsvt_data_one] : cd_.qsvt)
        {
            if(YMIX::compare_strings(qsvt_data_one.type, YVSv{"QSVT-ham", "QSP-ham"}))
            {
                counters_t[name_qsvt_gate] = -1;
            }
        }

        map<string, int> counters_gadget;
        for(auto const& [name_gadget, gadget_data] : cd_.gadgets)
        {
            if(YMIX::compare_strings(gadget_data.type, YVSv{"compression"}))
            {
                counters_gadget[name_gadget] = -1;
            }
        }

        timer_comp.Start();
        YMIX::print_log( "Calculating the circuit... ", 0, false, false);
        // printf("SIZE_MAX = %zu\n", SIZE_MAX);
        while(id_current_gate < u_work->get_N_gates())
        {
            // generate the circuit:
            // timer_comp.Start();
            // YMIX::print_log( "Calculating the circuit... ", 0, false, false);
            u_work->generate(stop_point_name, id_current_gate);

            if(YMIX::compare_strings(sel_compute_output_, "zero-ancillae"))
                u_work->get_state(outZ, true, false);
            if(YMIX::compare_strings(sel_compute_output_, "all"))
                u_work->get_state(outF);
            if(flag_zero_state_of_)
            {
                outZ_chosen.state_to_choose = qs_chosen_state_;
                u_work->get_state(outZ_chosen, false, false);
            }

            // timer_comp.Stop();
            // YMIX::print_log( "duration: " + timer_comp.get_dur_str_s());

            if(YMIX::compare_strings(sel_print_output_, "all"))
                YMIX::print_log(
                        "...Output all states after " + stop_point_name + ": \n" + outF.str_wv
                );
            else if(YMIX::compare_strings(sel_print_output_, "zero-ancillae"))
                YMIX::print_log(
                        "...Output zero-ancilla states after " + stop_point_name + ": \n" + outZ.str_wv
                );
            if(flag_zero_state_of_)
                YMIX::print_log(
                        "...Output chosen states after " + stop_point_name + ": \n" + outZ_chosen.str_wv
                );

            // --- Save results from each time step of QSP simulations ---
            if(stop_point_name.find("QSP-H") != string::npos && flag_output_qsp_)
            {
                auto pos1 = stop_point_name.find("<");
                auto pos2 = stop_point_name.find(">");
                string name_qsvt_gate = stop_point_name.substr(pos1+1, pos2-pos1-1);

                counters_t[name_qsvt_gate] += 1;

                u_work->get_state(outZ_qsp, true, true);

                hfo_.open_w();
                if(outZ_qsp.ampls.size() > 0)
                {
                    hfo_.add_vector(
                        outZ_qsp.ampls,  
                        "ampls-"s + to_string(counters_t[name_qsvt_gate]) + "-" + to_string(count_init_state), 
                        name_qsvt_gate
                    );
                }
                hfo_.close(); 
            }

            // --- Save results from GADGETS ---
            if(stop_point_name.find("CompressionGadget") != string::npos && flag_output_gadget_)
            {
                auto pos1 = stop_point_name.find("<");
                auto pos2 = stop_point_name.find(">");
                string name_gadget = stop_point_name.substr(pos1+1, pos2-pos1-1);

                counters_gadget[name_gadget] += 1;
                auto& gadget_data = cd_.gadgets[name_gadget];
                if(YMIX::compare_strings(gadget_data.type, "compression"))
                {
                    // --- Chose only states that contain correct state in the counter register ---
                    int des_index = gadget_data.N_mult - counters_gadget[name_gadget] - 1;
                    auto nc = gadget_data.counter_qubits.size();
                    auto n_full = u_work->get_n_qubits();
                    std::vector<short> des_bitstring(nc);
                    YMATH::intToBinary(des_index, des_bitstring);

                    std::vector<short> chosen_bitstring(n_full);
                    std::fill(chosen_bitstring.begin(), chosen_bitstring.end(), -1);
                    for(int iq = 0; iq < des_bitstring.size(); iq++)
                    {
                        auto one_bit = des_bitstring[nc-iq-1];
                        auto id_bit = gadget_data.counter_qubits[iq];
                        chosen_bitstring[n_full-id_bit-1] = one_bit;
                    }

                    if(flag_zero_state_of_)
                    {
                        for(int iq = 0; iq < chosen_bitstring.size(); iq++)
                            if(chosen_bitstring[iq] == -1)
                                chosen_bitstring[iq] = qs_chosen_state_[iq];
                        outZ_gadget.state_to_choose = chosen_bitstring;
                        u_work->get_state(outZ_gadget, false, false);
                        // cout << "FLAG: 1" << endl;
                    }
                    else
                    {
                        outZ_gadget.state_to_choose = chosen_bitstring;
                        u_work->get_state(outZ_gadget, false, false);
                        // cout << "FLAG: 2" << endl;
                    }
                    
                    hfo_.open_w();
                    if(outZ_gadget.ampls.size() > 0)
                    {
                        hfo_.add_vector(
                            outZ_gadget.ampls,  
                            "ampls-"s + to_string(counters_gadget[name_gadget]) + "-" + to_string(count_init_state), 
                            name_gadget
                        );
                        hfo_.add_matrix(
                            outZ_gadget.states, 
                            "states-"s + to_string(counters_gadget[name_gadget]) + "-" + to_string(count_init_state),     
                            name_gadget
                        );
                    }
                    hfo_.close();
                }
            }
        }

        // --- Store the output state at the very end of the circuit ---
        hfo_.open_w();
        if(YMIX::compare_strings(sel_compute_output_, "all"))
        {
            hfo_.add_vector(
                outF.ampls,  
                "output-all-amplitudes-"s + to_string(count_init_state), 
                "states"
            );
            hfo_.add_matrix(
                outF.states, 
                "output-all-states-"s + to_string(count_init_state),     
                "states"
            );
        }
        if(YMIX::compare_strings(sel_compute_output_, YVSv{"zero-ancillae"}))
            if(outZ.ampls.size() > 0)
            {
                hfo_.add_vector(
                    outZ.ampls,  
                    "output-zero-anc-amplitudes-"s + to_string(count_init_state), 
                    "states"
                );
                hfo_.add_matrix(
                    outZ.states, 
                    "output-zero-anc-states-"s + to_string(count_init_state),     
                    "states"
                );
            }
        if(flag_zero_state_of_)
            if(outZ_chosen.ampls.size() > 0)
            {
                hfo_.add_vector(
                    outZ_chosen.ampls,  
                    "chosen-amplitudes-"s + to_string(count_init_state), 
                    "states"
                );
                hfo_.add_matrix(
                    outZ_chosen.states, 
                    "chosen-states-"s + to_string(count_init_state),     
                    "states"
                );
                hfo_.add_vector(qs_chosen_state_,  "chosen_qubits"s, "states");
            }
        hfo_.close(); 
    }

    // --- calculate state probabilities on the indicated qubits ---
    {
        if(YMIX::compare_strings(sel_compute_output_, "none"))
        {
            timer_comp.Start();
            YMIX::print_log( "Calculating the circuit... ", 0, false, false);
            u_work->generate();
            timer_comp.Stop();
            YMIX::print_log( "duration: " + timer_comp.get_dur_str_s());
        }// otherwise, the circuit has already been generated;


        hfo_.open_w();
        if(flag_prob_)
        {
            vector<qreal> outProbs(1<<focus_qubits_.size());
            YMIX::print_log( "Calculating the probabilities... \n", 0, false, false);
            calcProbOfAllOutcomes(
                &outProbs[0], 
                u_work->get_qureg(), 
                &focus_qubits_[0],
                focus_qubits_.size()
            );

            hfo_.add_group("probabilities-"s + to_string(count_init_state));
            hfo_.add_vector(focus_qubits_,  "qubits"s, "probabilities-"s + to_string(count_init_state));
            hfo_.add_vector(outProbs,  "probs", "probabilities-"s + to_string(count_init_state));
        }
        hfo_.close();
    }

    timer_comp.Stop();
    YMIX::print_log( "duration: " + timer_comp.get_dur_str_s());
}


void QuCF__::calc_matrix(shared_ptr<QCircuit>& u_work, bool flag_sparse_format)
{
    // throughout the function, we assume that all nonancilla qubits are less significant
    // than any ancilla qubit;
    YMIX::print_log("\nCompute the circuit matrix");

    uint32_t n_all_qubits = u_work->get_n_qubits();
    uint32_t n_anc_qubits = u_work->get_na();
    uint32_t n_nonanc_qubits = n_all_qubits - n_anc_qubits;
    uint32_t N_matrix = 1 << n_nonanc_qubits;
    uint32_t N2 = N_matrix * N_matrix;

    // first register name -> register of a higher priority:
    vector<string> reg_names_nonanc;
    vector<int> n_qubits;
    u_work->get_nonancilla_regs(reg_names_nonanc, n_qubits);
    int N_reg = reg_names_nonanc.size();

    // take all possible input states in the nonancilla qubits and
    // compute the output states:
    YMIX::YTimer timer_comp;
    timer_comp.Start();

    YMATH::YMatrix A_real(N_matrix, N_matrix, !flag_sparse_format);
    YMATH::YMatrix A_imag(N_matrix, N_matrix, !flag_sparse_format);

    // one input state for each row index:
    for(uint32_t ir = 0; ir < N_matrix; ir++)
    {
        // reset the circuit:
        u_work->empty_binary_states();
        u_work->reset_qureg();

        // set an initial quantum state:
        vector<short> binArray(n_nonanc_qubits);
        YMATH::intToBinary(ir, binArray);

        // cout << "\n--------- row = " << ir << " ------------" << endl;

        int shift = 0;
        for(int ireg = 0; ireg < N_reg; ireg++)
        {
            string reg_name = reg_names_nonanc[ireg];
            int nq_reg      = n_qubits[ireg];
            auto reg_chosen = u_work->get_regs()[reg_name];

            // cout << "\n" << reg_name << " with " << nq_reg << endl;
            // cout << "chosen qubits: ";
            // for(auto ii = 0; ii < reg_chosen.size(); ii++)
            //     cout << reg_chosen[ii] << " ";

            vector<short> reg_bitstring(nq_reg);
            copy(
                binArray.begin() + shift, 
                binArray.begin() + shift + nq_reg,
                reg_bitstring.begin()
            );

            // cout << "bistring: ";
            // for(auto ii = 0; ii < reg_bitstring.size(); ii++)
            //     cout << reg_bitstring[ii] << " ";

            vector<int> reg_state;
            for(int id_bit = 0; id_bit < nq_reg; id_bit++)
                if(reg_bitstring[nq_reg - id_bit - 1] == 1)
                    reg_state.push_back(id_bit);

            // cout << "\nqubits: ";
            // for(auto ii = 0; ii < reg_state.size(); ii++)
            //     cout << reg_state[ii] << " ";
            // cout << endl;

            u_work->set_reg_state(reg_name, reg_state);
            shift += nq_reg;
        }
        u_work->set_init_binary_state();

        // compute an output state:
        // string stop_point_name;
        // int id_current_gate = 0;
        // u_work->generate(stop_point_name, id_current_gate);
        u_work->generate();

        // get the output state:
        YMIX::StateVectorOut out_state;
        u_work->get_state(out_state, true);

        // convert bitstrings to columns, 
        // convert amplitudes to values of matrix elements:
        int i_state = -1;
        for(auto one_state: out_state.states)
        {
            i_state++;
            Complex one_ampl = out_state.ampls[i_state];
            int id_column = YMATH::binaryToInt(one_state);

            // cout << "\n ---" << endl;
            // cout << "ic: " << id_column << endl;
            // cout << "ampl: " << one_ampl.real << " + " <<  one_ampl.imag << "*i" << endl;

            A_real(ir, id_column) = one_ampl.real;
            A_imag(ir, id_column) = one_ampl.imag;
        }
    }
    timer_comp.Stop();
    YMIX::print_log("duration: " + timer_comp.get_dur_str_s());
    YMIX::print_log("");

    // save the matrix to the .hdf5 file:
    hfo_.open_w();

    hfo_.add_group("matrix");
    hfo_.add_scalar(u_work->get_name(), "name-oracle", "matrix");
    hfo_.add_scalar(N_matrix, "N"s, "matrix"s);

    hfo_.add_array(A_real.get_1d_pointer(), N2, "real", "matrix"s);
    hfo_.add_array(A_imag.get_1d_pointer(), N2, "imag", "matrix"s);

    hfo_.close();
}