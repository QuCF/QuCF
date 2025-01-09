#include "../include/QCircuit.h"

using namespace std;




QCircuit::QCircuit(
    YCS name, 
    const QuESTEnv& env, 
    YCS path_to_output, 
    YCU nq, 
    const std::map<std::string, qreal>& constants,
    YCB flag_circuit,
    YCB flag_tex,
    YCB flag_layers,
    YCB flag_stop_gates,
    YCB flag_repeat_insert,
    YCB flag_progress_bar
) :
name_(name),
env_(env),
path_to_output_(path_to_output),
flag_circuit_(flag_circuit),
flag_tex_(flag_tex),
flag_layers_(flag_layers),
flag_stop_gates_(flag_stop_gates),
flag_repeat_insert_(flag_repeat_insert),
flag_progress_bar_(flag_progress_bar)
{
    timer_.Start();
    nq_ = 0;
    if(nq > 0)
        create(nq); 
    else   
        nq_ = nq; 

    constants_ = map<string, qreal>(constants);
}


void QCircuit::create(YCU nq)
{
    if(nq <= 0)
    {
        cerr << "Error: circuit " << name_ << ": number of qubits has to be greater than zero" << endl;
        exit(-1);
    }

    if(nq_ > 0)
    {
        cerr << "Error: circuit " << name_ << ": the circuit has been already created" << endl;
        exit(-1);
    }

    nq_ = nq;
    flag_circuit_allocated_ = false;
    // c_ = createQureg(nq_, env_);

    tex_noc_.resize(nq_, 1);
    tex_lines_.resize(nq_);
    tex_qubit_names_.resize(nq_);

    create_circ_file();
    create_tex_file();
    // initZeroState(c_);

    if(flag_layers_) oo_layers_ = make_shared<CircuitLayers__>(nq_);
}

void QCircuit::allocate_circuit()
{
    c_ = createQureg(nq_, env_);
    // initZeroState(c_);
    flag_circuit_allocated_ = true;
}


QCircuit::QCircuit(YCCQ oc, YCS cname)
{
    timer_.Start();

    if(cname.empty())
    {
        name_ = oc->name_ + "-copy";
        flag_circuit_ = false;
        flag_tex_     = false;
        flag_layers_  = false;
        flag_stop_gates_ = true;
        flag_repeat_insert_ = false;
        flag_progress_bar_ = false;
    }
    else
    {
        name_ = cname;
        flag_circuit_    = oc->flag_circuit_;
        flag_tex_        = oc->flag_tex_;
        flag_layers_     = oc->flag_layers_;
        flag_stop_gates_ = oc->flag_stop_gates_;
        flag_repeat_insert_ = oc->flag_repeat_insert_;
        flag_progress_bar_ = oc->flag_progress_bar_;
    }
        
    env_ = oc->env_;
    nq_ = oc->nq_;
    flag_circuit_allocated_ = false;
    // c_ = createQureg(nq_, env_);
    // initZeroState(c_);

    path_to_output_ = oc->path_to_output_;
    init_state_     = oc->init_state_;
    regs_           = map<string, YVIv>(oc->regs_);
    flags_anc_regs_ = map<string, bool>(oc->flags_anc_regs_);
    regnames_       = vector<string>(oc->regnames_);
    ancs_           = YVIv(oc->ancs_);
    ib_state_       = vector<short>(oc->ib_state_);
    id_start_       = oc->id_start_;
    standart_output_format_ = oc->standart_output_format_;
    // unique_gates_names_ = YVSv(oc->unique_gates_names_);

    // flag_circuit_ = oc->flag_circuit_;
    // flag_tex_     = oc->flag_tex_;
    // flag_layers_  = oc->flag_layers_;

    tex_noc_ = vector<uint64_t>(oc->tex_noc_);
    tex_lines_ = vector<vector<string>>(oc->tex_lines_);
    tex_qubit_names_ = vector<string>(oc->tex_qubit_names_);

    blocks_ids_unit_ = std::vector<YVIv>(oc->blocks_ids_unit_);
    blocks_ids_zero_       = std::vector<YVIv>(oc->blocks_ids_zero_);

    create_circ_file();
    create_tex_file();
    save_regs();

    for(auto const& gate: oc->gates_)
    {
        auto gate_copy = gate->copy_gate();
        gates_.push_back(gate_copy);
    }
    constants_ = map<string, qreal>(oc->constants_);
    if(flag_layers_) oo_layers_ = make_shared<CircuitLayers__>(oc->oo_layers_);
}


QCircuit::~QCircuit()
{
    finish_tex_file();
    if(flag_circuit_allocated_)
        destroyQureg(c_, env_);
    timer_.Stop();
    if(env_.rank == 0)
    {
        YMIX::LogFile cf;
        cf << "Circuit " << name_ << ": Life time is " << scientific << setprecision(3) << timer_.get_dur() << " ms\n";
    }
}


void QCircuit::create_circ_file()
{
    if(!flag_circuit_)
        return;

    cfname_ = path_to_output_ + "/" + name_ + FORMAT_CIRCUIT;
    if(env_.rank == 0)
    {
        YMIX::File cf(cfname_, true);
        cf << name_ << " " << nq_ << "\n";
    }
}


void QCircuit::create_tex_file()
{
    if(!flag_tex_)
    {
        // cout << "HERE: circuit [" << name_ << "]: no tex" << endl;
        return;
    }
        

    texname_ = path_to_output_ + "/" + name_ + FORMAT_TEX;
    YMIX::File cf(texname_, true);

    cf << "\\documentclass{article}\n";
    // cf << "\\usepackage[utf8]{inputenc}\n";
    cf << "\\usepackage[dvipsnames]{xcolor}\n";
    cf << "\\usepackage{amsmath}\n";
    cf << "\\usepackage{amsfonts,amssymb}\n";
    cf << "\\usepackage{tikz}\n";
    cf << "\\usetikzlibrary{quantikz}\n";
    cf << "\n";
    cf << "\\newcommand{\\yi}{\\mathrm{i}}\n";
    cf << "\n";
    cf << "\\begin{document}\n";
    cf << "\n\n";
    cf << "\\begin{figure}[t!]\n";
    cf << "\\centering\n";
    cf << TEX_BEGIN_CIRCUIT;   
}


void QCircuit::finish_tex_file()
{
    if(!flag_tex_)
        return;

    if(env_.rank == 0)
    {
        YMIX::File cf(texname_, false);

        // --- adjust the widths of the .tex phrases ---
        auto n_layers = tex_lines_[0].size();
        uint32_t width_max, width_curr;
        vector<uint32_t> q_widths(nq_);
        for(auto id_layer = 1; id_layer < n_layers; id_layer++)
        {
            if(id_layer%YGV::tex_circuit_length == 0)
            {
                for(auto id_q = 0; id_q < nq_; id_q++)
                    q_widths[id_q] = 0;
                continue;
            }

            for(auto id_q = 0; id_q < nq_; id_q++)
                q_widths[id_q] += tex_lines_[id_q][id_layer-1].length();

            width_max = 0;
            for(auto id_q = 0; id_q < nq_; id_q++)
                if(width_max < q_widths[id_q])
                    width_max = q_widths[id_q];

            for(auto id_q = 0; id_q < nq_; id_q++)
            {
                width_curr = q_widths[id_q];
                string line_padding(width_max - width_curr, ' ');

                tex_lines_[id_q][id_layer] = line_padding + tex_lines_[id_q][id_layer];
            }
        }

        // --- write the matrix of strings to the file ---
        int n_circ_pieces = int(tex_lines_[0].size() / YGV::tex_circuit_length)+1;
        if(tex_lines_[0].size() % YGV::tex_circuit_length == 0)
            n_circ_pieces--;
        for(auto id_piece = 0; id_piece < n_circ_pieces; id_piece++)
        {
            int counter_row = -1;
            int column_begin = id_piece * YGV::tex_circuit_length;
            int column_end = std::min(column_begin + YGV::tex_circuit_length, int(tex_lines_[0].size()));
            for(auto const & one_row_vec: tex_lines_)
            {
                counter_row++;
                string line_row = "";
                for(auto id_c = column_begin; id_c < column_end; id_c++)
                    line_row += one_row_vec[id_c];
                if(counter_row < (tex_lines_.size()-1))
                    line_row += "\\\\"s;
                
                if(id_piece > 0)
                    cf << tex_qubit_names_[counter_row];
                cf << line_row << "\n";
            }
            cf << "\\end{quantikz}";

            if(id_piece < (n_circ_pieces - 1))
            {
                cf << "\\\\\n";
                cf << TEX_VERTICAL_GAP << "\n";
                cf << TEX_BEGIN_CIRCUIT;
            }
            else
            {
                cf << "\n";
            }
        }

        // --- finish the .tex file ---
        string name_mod = name_;
        YMIX::replace_substrings(name_mod, "_", "\\_");

        cf << "\\caption{Figure of the circuit: " << name_mod << "}\n";
        cf << "\\end{figure}\n";
        cf << "\n\n";
        cf << "\\end{document}\n";
    }
}


void QCircuit::generate()
{
    if(!flag_circuit_allocated_)
        throw string("The circuit [" + name_ + "] is not allocated in memory.");

    YMIX::YTimer timer;
    auto start = gates_.begin() + id_start_;
    for(auto it = start; it != gates_.end(); ++it)
    {
        (*it)->generate(c_);
    }
    id_start_ += gates_.end() - start;
}


void QCircuit::generate(string& stop_name, uint64_t& id_current)
{
    if(!flag_circuit_allocated_)
        throw string("The circuit [" + name_ + "] is not allocated in memory.");

    auto start  = gates_.begin() + id_start_;
    auto it_end = gates_.end();
    stop_name = name_;

    int perc_prev = 0;
    float perc_float = 0.;
    int perc_int = 0;
    int Ng = gates_.size();
    int counter_g = -1;
    int counter_ratio = 1;
    int step_prog = 2;
    if(id_start_ == 0 & !flag_stop_gates_)
        cout << "\nProgress:\n" << std::flush;
    for(auto it = start; it != gates_.end(); ++it) 
    {
        if(!flag_stop_gates_)
        {
            counter_g++;
            perc_float = (counter_g*1.)/(Ng-1) * 100.;

            if(flag_progress_bar_)
            { 
                perc_int = perc_float;
                if((perc_int - perc_prev) >= step_prog)
                {
                    for(int ii = 0; ii < 12; ii++)
                        cout << "\b";
                    cout << "| " << std::setw(10) << perc_float << "%" << std::flush;
                    perc_prev = perc_int;
                }
                else
                {
                    if(counter_g%1000 == 0)
                    {
                        if(counter_g > 0)
                            for(int ii = 0; ii < 12; ii++)
                                cout << "\b";
                        cout << std::setw(11) << perc_float << "%" << std::flush;
                    }
                }
                if(perc_int == 100)
                    cout << endl;
            }
            else
            {   
                if(perc_float >= counter_ratio * 5.)
                {
                    printf("%6.1f...", perc_float);
                    if(counter_ratio%10 == 0)
                        cout << std::endl;
                    else
                        cout << std::flush;
                    counter_ratio++;
                }
            }
        }

        if(YMIX::compare_strings((*it)->get_type(), "stop"))
        {
            it_end = it+1;
            stop_name = (*it)->get_name();
            break;
        }

        (*it)->generate(c_);
    }
    id_start_ += it_end - start;
    id_current = id_start_;
}


void QCircuit::print_gates()
{
    // -------------------------------------------------
    // --- print gates to the .circuit file ---
    if(env_.rank == 0 && flag_circuit_)
    {
        YMIX::File cf(cfname_);
        for(auto& gate: gates_)
            gate->write_to_file(cf);
    }

    // -------------------------------------------------
    // --- print gates to the .tex file ---
    if(env_.rank == 0 && flag_tex_)
    {
        YVIv ids_qubits_of_gate, ids_range, ids_targets;
        uint64_t id_first_noc_layer;
        uint64_t id_new_noc_layer;
        int id_b, id_t;
        bool flag_box = false;
        string box_name;
        for(auto& gate: gates_)
        {
            id_first_noc_layer = 0;
            if(YMIX::compare_strings(gate->get_type(), "stop"))
                continue;

            if(YMIX::compare_strings(gate->get_type(), "box"))
            {
                string box_name_curr = gate->get_name();
                if(gate->get_flag_start())
                {
                    if(!flag_box)
                        box_name = box_name_curr;
                    flag_box = true;
                }
                else
                {
                    if(YMIX::compare_strings(box_name, box_name_curr))
                        flag_box = false;
                }  
            }
                
            if(flag_box)
                continue;

            gate->get_gubits_act_on(ids_qubits_of_gate);
            id_t = *(max_element(ids_qubits_of_gate.begin(), ids_qubits_of_gate.end()));
            id_b = *(min_element(ids_qubits_of_gate.begin(), ids_qubits_of_gate.end()));

            // --- find the first layer in .tex, where there is enough free qubits to place the gate ---
            for(auto id_qubit = id_b; id_qubit <= id_t; id_qubit++)
            {
                if(id_first_noc_layer < tex_noc_[nq_ - id_qubit - 1])
                    id_first_noc_layer = tex_noc_[nq_ - id_qubit - 1];
            }

            // --- put the gate to the .tex layer ---
            gate->write_tex(tex_lines_, id_first_noc_layer, nq_);

            // --- shift the ids of non-occupied layers in the .tex ---
            id_new_noc_layer = id_first_noc_layer + 1;
            for(auto id_qubit = id_b; id_qubit <= id_t; id_qubit++)
                tex_noc_[nq_ - id_qubit - 1] = id_new_noc_layer;

            // --- if necessary, add next empty .tex layer ---
            if(id_new_noc_layer >= tex_lines_[0].size())
                for(auto id_q = 0; id_q < nq_; id_q++)
                    tex_lines_[id_q].push_back("&\\qw"s);
        }
    }
}


void QCircuit::h_adjoint()
{
    reverse(gates_.begin(), gates_.end());
    for(auto& gate: gates_)
    {
        gate->h_adjoint();
        if(flag_layers_) gate->set_layer(oo_layers_->get_n_layers() - gate->get_layer());
    }  
}


void QCircuit::copy_gates_from(
    YCCQ c, YCVI regs_new, YCCB box, YCVI cs_unit, YCVI cs_zero, YCB flag_inv
){
    if(box)
    {
        YSG oo = box->copy_gate();
        oo->add_control_qubits(cs_unit, cs_zero);
        if(flag_layers_) oo_layers_->add_gate(oo);
        gates_.push_back(oo);
    }

    if(flag_inv)
    {
        // auto gates_c = vector<YSG>(c->gates_);
        // reverse(gates_c.begin(), gates_c.end());
        // for(auto& gate: gates_c)
        // {
        //     auto gate_copy = gate->copy_gate();

        //     gate_copy->h_adjoint();
        //     gate_copy->correct_qubits(regs_new);
        //     gate_copy->add_control_qubits(cs_unit, cs_zero);
        //     if(flag_layers_) oo_layers_->add_gate(gate_copy);
        //     gates_.push_back(gate_copy);
        // }

        for(auto gate_it = c->gates_.rbegin(); gate_it != c->gates_.rend(); ++gate_it)
        {
            auto gate_copy = (*gate_it)->copy_gate();
            gate_copy->h_adjoint();
            gate_copy->correct_qubits(regs_new);
            gate_copy->add_control_qubits(cs_unit, cs_zero);
            if(flag_layers_) oo_layers_->add_gate(gate_copy);
            gates_.push_back(gate_copy);
        }
    }
    else
    {
        for(auto& gate: c->gates_)
        {
            auto gate_copy = gate->copy_gate();
            gate_copy->correct_qubits(regs_new);
            gate_copy->add_control_qubits(cs_unit, cs_zero);
            if(flag_layers_) oo_layers_->add_gate(gate_copy);
            gates_.push_back(gate_copy);
        }
    }

    if(box)
    {
        YSG oo = box->copy_box();
        oo->set_flag_start(false);
        oo->add_control_qubits(cs_unit, cs_zero);
        if(flag_layers_) oo_layers_->add_gate(oo);
        gates_.push_back(oo);
    }

    // unique_gates_names_.insert(
    //     unique_gates_names_.end(), 
    //     c->unique_gates_names_.begin(), 
    //     c->unique_gates_names_.end()
    // );
}


void QCircuit::insert_gates_from(const QCircuit* c, YCCB box)
{
    if(box)
    {
        YSG oo = box->copy_gate();
        if(flag_layers_) oo_layers_->add_gate(oo);
        gates_.push_back(oo);
    }


    for(const auto& gate: c->gates_)
    {
        // !!! here, do not add the inserted gates to the layers !!!
        gates_.push_back(gate);
    }

        
    if(box)
    {
        YSG oo = box->copy_box();
        oo->set_flag_start(false);
        if(flag_layers_) oo_layers_->add_gate(oo);
        gates_.push_back(oo);
    }

    // unique_gates_names_.insert(
    //     unique_gates_names_.end(), 
    //     c->unique_gates_names_.begin(), 
    //     c->unique_gates_names_.end()
    // );
}


YVIv QCircuit::add_register(YCS name, YCU n_qubits, YCB flag_ancilla)
{
    unsigned shift = 0;
    for(auto const& reg_name: regnames_)
        shift += regs_[reg_name].size();

    regnames_.push_back(name);

    // 0-th qubit is the least significant in the register;
    // the last qubit is the most significant:
    YVIv qubits_positions(n_qubits);
    for(unsigned i = 0; i < n_qubits; ++i)
        qubits_positions[i] = nq_ - shift - n_qubits + i; // YINV

    regs_[name] = qubits_positions;

    if(flag_ancilla)
    {
        ancs_.insert(
            ancs_.begin(), 
            qubits_positions.begin(), 
            qubits_positions.end()
        );
        flags_anc_regs_[name] = true;
    }
    else
    {
        flags_anc_regs_[name] = false;
    }
        
    return qubits_positions;
}


void QCircuit::print_reg_positions(std::ofstream& of) const
{
    for(auto const& reg_name: regnames_)
    {
        of << reg_name << " ";
        for(auto& id_q: regs_.at(reg_name))
            of << id_q << " ";
    }
}


void QCircuit::set_standart_output_format()
{
    int reg_size;
    int reg_size_total = 0;
    if(regs_.empty())
    {
        standart_output_format_.push_back(nq_);
        return;
    }
    
    for(auto const& reg_name: regnames_)
    {
        reg_size = regs_[reg_name].size();
        reg_size_total += reg_size;
        standart_output_format_.push_back(reg_size);
    }

    int diff = nq_ - reg_size_total;
    if(diff > 0) standart_output_format_.push_back(diff);
}


void QCircuit::save_regs()
{
    if(env_.rank == 0 && flag_circuit_) 
    {
        if(regnames_.size() > 0)
        {
            YMIX::File cf(cfname_);
            cf << "QubitRegisterNames ";
            for(auto const& reg_name: regnames_)
            {
                cf << reg_name << " " << regs_[reg_name].size() << " ";
                if(flags_anc_regs_[reg_name])
                    cf << "1" << " ";
                else
                    cf << "0" << " ";
            }
            cf.of << endl;
        }
    }

    if(env_.rank == 0 && flag_tex_) 
    {
        int counter_q = -1;
        for(auto const& reg_name: regnames_)
        {
            auto reg_qs = regs_[reg_name];
            auto reg_nq = reg_qs.size();

            string reg_name_part1, reg_name_part2;
            std::size_t pos = reg_name.find("_");
            if(pos != std::string::npos)
            {
                reg_name_part1 = reg_name.substr(0, pos);
                reg_name_part2 = reg_name.substr(pos+1) + ", ";
            }
            else
            {
                reg_name_part1 = string(reg_name);
                reg_name_part2 = "";
            }

            for(auto i = 0; i < reg_nq; i++)
            {
                counter_q++;

                string qu_name = "\\lstick{$" + reg_name_part1 + "_{" + reg_name_part2 + to_string(reg_nq - i - 1) + "}$}";
                tex_qubit_names_[counter_q] = qu_name;

                // register name and qubit id within the register:
                tex_lines_[counter_q].push_back(qu_name);

                // first layer:
                tex_lines_[counter_q].push_back("&\\qw");
            }
        }
    }
}


YQCP QCircuit::get_the_circuit()
{
    YQCP temp = this;
    return this;
}

void QCircuit::empty_binary_states()
{
    ib_state_ = vector<short>(nq_);
}
void QCircuit::prepare_zero_init_state()
{
    if(!flag_circuit_allocated_)
        throw string("The circuit [" + name_ + "] is not allocated in memory.");
    initZeroState(c_);
}
void QCircuit::reset()
{
    gates_.clear();
    id_start_ = 0;

    blocks_ids_unit_.clear();
    blocks_ids_zero_.clear();

    if(flag_circuit_allocated_)
    {
        destroyQureg(c_, env_);
        c_ = createQureg(nq_, env_);
        if(init_state_.flag_defined)
            reset_init_vector(init_state_);
        else{
            if(!ib_state_.empty())
                set_init_binary_state();
        }
    }
    // unique_gates_names_.clear();
}
void QCircuit::reset_qureg()
{
    id_start_ = 0;
    if(flag_circuit_allocated_)
    {
        destroyQureg(c_, env_);
        allocate_circuit();
    } 
    else
        allocate_circuit();
}

void QCircuit::set_init_binary_state(const bool& flag_mpi_bcast)
{
    if(!flag_circuit_allocated_)
        throw string("The circuit [" + name_ + "] is not allocated in memory.");
    if(flag_mpi_bcast)
    {
        if(env_.rank > 0)
            if(ib_state_.empty())
                ib_state_ = vector<short>(nq_);
    }
    long long int ii = YMATH::binaryToInt(ib_state_);
    initClassicalState(c_, ii);
}


void QCircuit::set_init_vector(YVQ ampl_vec_real, YVQ ampl_vec_imag)
{
    long long b_ampl = 0;
    long long N_ampls = ampl_vec_real.size();
    if(N_ampls != ampl_vec_imag.size())
        throw string("vectors with real and imaginary parts of the initial state are of different size.");
    setAmps(c_, b_ampl, &ampl_vec_real[0], &ampl_vec_imag[0], N_ampls);

    init_state_.flag_defined = true;
    init_state_.b_ampl = b_ampl;
    init_state_.n_ampls = N_ampls;
    init_state_.ampl_vec_real = YVQv(ampl_vec_real);
    init_state_.ampl_vec_imag = YVQv(ampl_vec_imag);
}


void QCircuit::reset_init_vector(INIT_STATE__& state)
{
    setAmps(
        c_, state.b_ampl, 
        &state.ampl_vec_real[0], &state.ampl_vec_imag[0], 
        state.n_ampls
    );
}


void QCircuit::set_qubit_state(YCU id_q)
{
    if(ib_state_.empty())
        ib_state_ = vector<short>(nq_);
        
    if(id_q >= nq_)
    {
        cerr << "Error while setting init. state of qubits:\n";
        cerr << "A qubit with id = " << id_q << " is requested while there are only " 
            << nq_ << " qubits in a circuit " << name_ << "\n";
        cerr << "Remark: qubit indices start from zero." << endl;
        exit(-1);
    }
    ib_state_[nq_ - id_q - 1] = 1;
}

void QCircuit::set_qubit_state(YCVI ids_qs)
{
    if(ib_state_.empty())
        ib_state_ = vector<short>(nq_);

    for(auto const& id_q: ids_qs)
        set_qubit_state(id_q);
}

void QCircuit::set_reg_state(YCS name, YCI id_reg_qubit)
{
    if(ib_state_.empty())
        ib_state_ = vector<short>(nq_);

    if(YMIX::compare_strings(name, YGV::reg_whole_circuit))
    {
        set_qubit_state(id_reg_qubit);
        return;
    }

    if(regs_.find(name) == regs_.end())
    {
        cerr << "\nError: No register with the name " << name << " in a circuit " << name_ << "." << endl;
        exit(-1);
    }

    auto reg_qubits = regs_[name];
    if(id_reg_qubit >= reg_qubits.size())
    {
        cerr << "\nError: The register (" << name << ") has only "
            << reg_qubits.size() << " qubits,\n";
        cerr << "while a qubit with id = " << id_reg_qubit  
            << " has been requested." << endl;
        exit(-1);
    }

    int id_circuit_qubit = reg_qubits[id_reg_qubit];
    ib_state_[nq_ - id_circuit_qubit - 1] = 1;
}

void QCircuit::set_reg_state(YCS name, YCVI ids_reg_qubits)
{
    if(YMIX::compare_strings(name, YGV::reg_whole_circuit))
    {
        set_qubit_state(ids_reg_qubits);
        return;
    }

    if(ib_state_.empty())
        ib_state_ = vector<short>(nq_);

    for(auto& id_reg_qubit: ids_reg_qubits)
        set_reg_state(name, id_reg_qubit);
}

YQCP QCircuit::x(YCVI ts, YCVI cs_unit, YCVI cs_zero)
{ 
    for(auto& t:ts) x(t, cs_unit, cs_zero);
    return get_the_circuit();
}
YQCP QCircuit::y(YCVI ts, YCVI cs_unit, YCVI cs_zero)
{ 
    for(auto& t:ts) y(t, cs_unit, cs_zero);
    return get_the_circuit();
}
YQCP QCircuit::z(YCVI ts, YCVI cs_unit, YCVI cs_zero)
{ 
    for(auto& t:ts) z(t, cs_unit, cs_zero);
    return get_the_circuit();
}
YQCP QCircuit::h(YCVI ts, YCVI cs_unit, YCVI cs_zero)
{
    for(auto& t:ts) h(t, cs_unit, cs_zero);
    return get_the_circuit();
}


void QCircuit::get_ref_to_state_vector(qreal*& state_real, qreal*& state_imag)
{
    if(!flag_circuit_allocated_)
        throw string("The circuit [" + name_ + "] is not allocated in memory.");
    copyStateFromGPU(c_);
    state_real = c_.stateVec.real;
    state_imag = c_.stateVec.imag;
}


void QCircuit::get_state(
    YMIX::StateVectorOut& out, 
    YCB flag_ZeroPriorAnc, 
    YCB flag_zero_ampls 
){
    if(!flag_circuit_allocated_)
        throw string("The circuit [" + name_ + "] is not allocated in memory.");
    out.n_low_prior_qubits = flag_ZeroPriorAnc ? (nq_ - ancs_.size()): nq_;
    out.organize_state = get_standart_output_format();
    YMIX::Wavefunction_Probabilities(c_, out, flag_zero_ampls);
}


void QCircuit::read_structure_gate(
    YISS istr, YVI ids_target, qreal& par_gate, YVI ids_unit, YVI ids_zero
){
    string word;

    // --- read target qubits ---
    read_reg_int(istr, ids_target);
    
    // --- read the parameter of the gate ---  
    try
    {
        if(!isnan(par_gate))
        {
            istr >> word;
            par_gate = get_value_from_word(word);
        } 
    }   
    catch(YCS e)
    {
        throw "error in the format of the gate parameter: oracletool sees [" + word + "]: " + e;
    }

    // --- read the end of the gate structure description --- 
    read_end_gate(istr, ids_unit, ids_zero);
}


void QCircuit::read_structure_gate(
    YISS istr, YVI ids_target, qreal& par_gate1, qreal& par_gate2,
    YVI ids_unit, YVI ids_zero
){
    string word;

    // --- read target qubits ---
    read_reg_int(istr, ids_target);
    
    // --- read parameters of the gate ---  
    try
    {   
        if(!isnan(par_gate1))
        {
            istr >> word;
            par_gate1 = get_value_from_word(word);
        } 
        if(!isnan(par_gate2))
        {
            istr >> word;
            par_gate2 = get_value_from_word(word);
        } 
    }
    catch(YCS e)
    {
        throw "error in the format of the gate parameter: oracletool sees [" + word + "]: " + e;
    }

    // --- read the end of the gate structure description ---
    read_end_gate(istr, ids_unit, ids_zero);  
}


void QCircuit::read_end_element(YISS istr, YVI ids_unit, YVI ids_zero, YCU id_element)
{
    std::string word;

    if(id_element != 2)
    {
        for(YVIv& ids_zero_global: blocks_ids_zero_)
            ids_zero.insert(
                ids_zero.end(), 
                ids_zero_global.begin(), 
                ids_zero_global.end()
            );
        for(YVIv& ids_unit_global: blocks_ids_unit_)
            ids_unit.insert(
                ids_unit.end(), 
                ids_unit_global.begin(), 
                ids_unit_global.end()
            );
    }
    
    while(istr >> word)
    {
        if(id_element == 0)
            if(YMIX::compare_strings(word, "end_gate"))
                return;
        if(id_element == 1)
            if(YMIX::compare_strings(word, "end_circuit"))
            {
                // cout << "heree" << endl;
                return;
            }
        if(id_element == 2)
            if(YMIX::compare_strings(word, "do"))
                return;

        if(YMIX::compare_strings(word, "control"))
            read_reg_int(istr, ids_unit);
        if(YMIX::compare_strings(word, "ocontrol"))
            read_reg_int(istr, ids_zero);
        if(YMIX::compare_strings(word, "control_e"))
            read_reg_int(istr, ids_unit, ids_zero);
        if(YMIX::compare_strings(word, "ocontrol_e"))
            read_reg_int(istr, ids_zero, ids_unit);
    }
}


void QCircuit::read_reg_int(YISS istr, YVI ids_target, YCB flag_sort, YCS word_start)
{
    YVIv ids_target_e;
    read_reg_int_CORE(istr, ids_target, ids_target_e, flag_sort, word_start, false);
}


void QCircuit::read_reg_int(YISS istr, YVI ids_target, YVI ids_target_e, YCB flag_sort, YCS word_start)
{
    read_reg_int_CORE(istr, ids_target, ids_target_e, flag_sort, word_start, true);
}


void QCircuit::read_reg_int_CORE(
    YISS istr, YVI ids_target, YVI ids_target_e, 
    YCB flag_sort, YCS word_start, YCB flag_e
){
    string reg_name, word;
    bool flag_read_reg_name = false;
    int n_regs, integer_qu, n_bitA;
    int nq_reg;
    size_t pos1 = string::npos;

    if(word_start.empty())
        istr >> word;
    else 
        word = word_start;

    if(YMATH::is_number(word))
    {
        n_regs = stoi(word);
        flag_read_reg_name = true;
    }
    else
    {
        n_regs = 1;
        reg_name = word;
        flag_read_reg_name = false;
    }

    // within every register, one can have several qubits
    for(unsigned i_reg = 0; i_reg < n_regs; i_reg++)
    {
        if(flag_read_reg_name) 
            istr >> reg_name;

        pos1 = reg_name.find("[",0);
        if(pos1 != string::npos)
        {
            word = reg_name.substr(pos1);
            reg_name = reg_name.substr(0, size(reg_name) - size(word));
        }
        else
            word = "";
   
        if(!YMIX::is_present(regnames_, reg_name))
            throw "no register with the name " + reg_name;

        try
        {
            bool flag_empty = false;

            auto reg_chosen = regs_[reg_name];
            nq_reg = reg_chosen.size();
            if(word.empty())
                istr >> word;
            pos1 = word.find("[",0);

            // read an integer -> convert into a bistring 
            // -> target qubits are qubits in the unit states;
            if(pos1 == string::npos)
            {
                integer_qu = get_value_from_word(word);
                if(integer_qu < 0)
                    integer_qu = (1 << nq_reg) + integer_qu;

                vector<short> binArray(nq_reg);
                YMATH::intToBinary(integer_qu, binArray);

                n_bitA = binArray.size();
                for(unsigned id_bit = 0; id_bit < n_bitA; id_bit++)
                    if(binArray[n_bitA - id_bit - 1] == 1)
                        ids_target.push_back(reg_chosen[id_bit]);
                    else if(flag_e)
                        ids_target_e.push_back(reg_chosen[id_bit]);
            }
            // read an array of qubits;
            // the least significant qubit in a register has an index = 0;
            else
            {
                word = word.substr(pos1+1); // remove "["
                if(word.empty())
                    istr >> word;

                auto pos2 = word.find("]",0);
                int id_qu;
                vector<int> array_qu_e = vector<int>(reg_chosen);
                while(pos2 == string::npos)
                {
                    get_id_qu_pattern(id_qu, word, nq_reg);
                    ids_target.push_back(reg_chosen[id_qu]);
                    if(flag_e)
                        array_qu_e.erase(
                            find(array_qu_e.begin(), array_qu_e.end(), reg_chosen[id_qu])
                        );

                    istr >> word;
                    pos2 = word.find("]",0);
                }
                word = word.substr(0, size(word)-1); // remove "]"
                if(word.empty())
                    flag_empty = true;

                if(!flag_empty)
                {
                    get_id_qu_pattern(id_qu, word, nq_reg);
                    ids_target.push_back(reg_chosen[id_qu]);
                    if(flag_e) 
                        array_qu_e.erase(
                            find(array_qu_e.begin(), array_qu_e.end(), reg_chosen[id_qu])
                        );
                }
                
                if(flag_e)
                {
                    for(auto qu_e: array_qu_e)
                        ids_target_e.push_back(qu_e);
                }
            }
        }
        catch(YCS e)
        {
            throw "wrong format of the number of qubits in the register " + reg_name + ": " + e;
        }
    }

    // sort the final array:
    if(flag_sort)
        sort(ids_target.begin(), ids_target.end());
}


void QCircuit::read_structure_gate_adder_subtractor(
    YISS istr, YCS path_in, YCB flag_inv, YCI gate_type
){
    YVIv ids_target, ids_unit, ids_zero;
    long long nt;

    // --- read target qubits ---
    read_reg_int(istr, ids_target);

    // --- read end of the gate structure ---
    read_end_gate(istr, ids_unit, ids_zero);

    // add the operator:
    if(gate_type == 1) adder_1(ids_target, ids_unit, ids_zero, flag_inv);
    if(gate_type == 2) adder_2(ids_target, ids_unit, ids_zero, flag_inv);
    if(gate_type == 3) adder_3(ids_target, ids_unit, ids_zero, flag_inv);

    if(gate_type == -1) subtractor_1(ids_target, ids_unit, ids_zero, flag_inv);
    if(gate_type == -2) subtractor_2(ids_target, ids_unit, ids_zero, flag_inv);
    if(gate_type == -3) subtractor_3(ids_target, ids_unit, ids_zero, flag_inv);
}


void QCircuit::read_structure_gate_adder(YISS istr, YCS path_in, YCB flag_inv)
{
    YVIv ids_target_1, ids_target_2, ids_target_3, ids_unit, ids_zero;
    long long nt;

    // --- read target qubits ---
    read_reg_int(istr, ids_target_1);
    read_reg_int(istr, ids_target_2);
    read_reg_int(istr, ids_target_3);

    if(ids_target_1.size() != ids_target_2.size())
        throw string("target registers must have the same size");
    if(ids_target_1.size() != ids_target_3.size())
        throw string("target registers must have the same size");
    if(ids_target_2.size() != ids_target_3.size())
        throw string("target registers must have the same size");

    // --- read the end of the gate structure ---
    read_end_gate(istr, ids_unit, ids_zero);

    // --- adder of two variables ---
    adder(ids_target_1, ids_target_2, ids_target_3, ids_unit, ids_zero, flag_inv);
}


void QCircuit::read_structure_gate_subtractor(YISS istr, YCS path_in, YCB flag_inv)
{
    YVIv ids_target_1, ids_target_2, ids_target_3, ids_unit, ids_zero;
    long long nt;

    // --- read target qubits ---
    read_reg_int(istr, ids_target_1);
    read_reg_int(istr, ids_target_2);
    read_reg_int(istr, ids_target_3);

    if(ids_target_1.size() != ids_target_2.size())
        throw string("target registers must have the same size");
    if(ids_target_1.size() != ids_target_3.size())
        throw string("target registers must have the same size");
    if(ids_target_2.size() != ids_target_3.size())
        throw string("target registers must have the same size");

    // --- read the end of the gate structure ---
    read_end_gate(istr, ids_unit, ids_zero);

    // --- adder of two variables --- 
    subtractor(ids_target_1, ids_target_2, ids_target_3, ids_unit, ids_zero, flag_inv);
}


void QCircuit::read_structure_gate_adder_qft(YISS istr, YCS path_in, YCB flag_inv)
{
    YVIv qs_v1, qs_v2, qs_carry, ids_unit, ids_zero;
    int q_carry;
    long long nt;

    // --- read target qubits ---
    read_reg_int(istr, qs_v1);
    read_reg_int(istr, qs_v2);
    read_reg_int(istr, qs_carry);
    q_carry = qs_carry[0];

    if(qs_v1.size() != qs_v2.size())
        throw string("target registers must have the same size");

    // --- read the end of the gate structure ---
    read_end_gate(istr, ids_unit, ids_zero);

    // --- adder ---
    adder_qft(qs_v1, qs_v2, q_carry, ids_unit, ids_zero, flag_inv);
}


void QCircuit::read_structure_gate_subtractor_qft(YISS istr, YCS path_in, YCB flag_inv)
{
    YVIv qs_v1, qs_v2, qs_sign, ids_unit, ids_zero;
    int q_sign;
    long long nt;

    // --- read target qubits ---
    read_reg_int(istr, qs_v1);
    read_reg_int(istr, qs_v2);
    read_reg_int(istr, qs_sign);
    q_sign = qs_sign[0];

    if(qs_v1.size() != qs_v2.size())
        throw string("target registers must have the same size");

    // --- read the end of the gate structure ---
    read_end_gate(istr, ids_unit, ids_zero);

    // --- subtractor ---
    subtractor_qft(qs_v1, qs_v2, q_sign, ids_unit, ids_zero, flag_inv);
}



void QCircuit::read_structure_gate_adder_fixed(YISS istr, YCS path_in, YCB flag_inv)
{
    YVIv ids_target, ids_carry_temp, ids_unit, ids_zero;
    int  id_carry;
    uint32_t int_sub;
    string word;

    // --- read the gate parameters ---
    read_reg_int(istr, ids_target);
    
    istr >> word;
    int_sub = get_value_from_word(word);

    read_reg_int(istr, ids_carry_temp);
    id_carry = ids_carry_temp[0];

    // --- read the end of the gate structure ---
    read_end_gate(istr, ids_unit, ids_zero);

    // --- adder of two variables ---
    adder_fixed(ids_target, id_carry, int_sub, ids_unit, ids_zero, flag_inv);
}


void QCircuit::read_structure_gate_subtractor_fixed(YISS istr, YCS path_in, YCB flag_inv)
{
    YVIv ids_target, ids_carry_temp, ids_unit, ids_zero;
    int  id_carry;
    uint32_t int_sub;
    string word;

    // --- read the gate parameters ---
    read_reg_int(istr, ids_target);
    
    istr >> word;
    int_sub = get_value_from_word(word);

    read_reg_int(istr, ids_carry_temp);
    id_carry = ids_carry_temp[0];

    // --- read the end of the gate structure ---
    read_end_gate(istr, ids_unit, ids_zero);

    // --- gate---
    subtractor_fixed(ids_target, id_carry, int_sub, ids_unit, ids_zero, flag_inv);
}


void QCircuit::read_structure_gate_comparator_fixed(YISS istr, YCS path_in, YCB flag_inv)
{
    YVIv ids_target, ids_carry, ids_unit, ids_zero;
    uint32_t int_sub;
    string word;

    // --- read the gate parameters ---
    read_reg_int(istr, ids_target);
    
    istr >> word;
    int_sub = get_value_from_word(word);

    read_reg_int(istr, ids_carry);

    // --- read the end of the gate structure ---
    read_end_gate(istr, ids_unit, ids_zero);

    // --- gate ---
    comparator_fixed(ids_target, ids_carry, int_sub, ids_unit, ids_zero, flag_inv);
}


void QCircuit::read_structure_gate_swap(YISS istr, YCS path_in, YCB flag_inv)
{
    YVIv ids_target_1, ids_target_2, ids_unit, ids_zero;
    long long nt;

    // --- read target qubits ---
    read_reg_int(istr, ids_target_1);
    read_reg_int(istr, ids_target_2);

    // --- read end of gate structure ---
    read_end_gate(istr, ids_unit, ids_zero);

    // --- add CNOT gates ---
    nt = ids_target_1.size();
    for(unsigned i = 0; i < nt; ++i)
        swap(ids_target_1[i], ids_target_2[i], ids_unit, ids_zero);
}


void QCircuit::read_structure_gate_fourier(YISS istr, YCS path_in, YCB flag_inv)
{
    YVIv ids_target, ids_unit, ids_zero;

    // --- read target qubits ---
    read_reg_int(istr, ids_target);

    // --- read end of gate structure ---
    read_end_gate(istr, ids_unit, ids_zero);

    // add the quantum Fourier circuit:
    quantum_fourier(ids_target, ids_unit, ids_zero, flag_inv, true);
}


void QCircuit::read_structure_sin(YISS istr, YCS path_in, YCB flag_inv)
{
    YVIv ids_a, ids_main, ids_unit, ids_zero;
    qreal alpha_0, alpha;
    string word;

    // --- read an ancilla qubit to put rotations there ---
    read_reg_int(istr, ids_a);

    // --- read angles ---
    istr >> word;
    alpha_0 = get_value_from_word(word);

    istr >> word;
    alpha = get_value_from_word(word);

    // --- read the condition qubits ---
    read_reg_int(istr, ids_main);

    // --- read end of gate structure ---
    read_end_gate(istr, ids_unit, ids_zero);

    // add the quantum Fourier circuit:
    gate_sin(ids_a, ids_main, alpha_0, alpha, ids_unit, ids_zero, flag_inv);
}


void QCircuit::read_structure_sinC(YISS istr, YCS path_in, YCB flag_inv)
{
    YVIv ids_a, ids_main, ids_unit, ids_zero;
    qreal alpha_0_y, alpha_y, alpha_0_z, alpha_z;
    string word;

    // --- read an ancilla qubit to put rotations there ---
    read_reg_int(istr, ids_a);

    // --- read angles ---
    istr >> word;
    alpha_0_y = get_value_from_word(word);

    istr >> word;
    alpha_y = get_value_from_word(word);

    istr >> word;
    alpha_0_z = get_value_from_word(word);

    istr >> word;
    alpha_z = get_value_from_word(word);

    // --- read the condition qubits ---
    read_reg_int(istr, ids_main);

    // --- read end of gate structure ---
    read_end_gate(istr, ids_unit, ids_zero);

    // add the quantum Fourier circuit:
    gate_sinC(
        ids_a, ids_main, 
        alpha_0_y, alpha_y, alpha_0_z, alpha_z, 
        ids_unit, ids_zero, flag_inv
    );
}


void QCircuit::read_structure_compression_gadget(
    YISS istr, 
    std::map<std::string, 
    YSQ>& ocs, 
    YCB flag_inv, 
    QuCF_complex_data& qucf_d
){
    string gadget_name;
    YVIv ids_counter;
    YVIv ids_U_target;
    string name_U;
    int N_mult;
    YVIv ids_unit, ids_zero; 
    bool flag_step_output;
    string word;

    auto& map_gadget_data = qucf_d.gadgets;
    GADGET_pars data;

    bool flag_data_read;

    // --- read the unique name of the gadget ---
    istr >> gadget_name;
    flag_data_read = qucf_d.check_name(gadget_name);

    data.name = gadget_name;
    data.type = "compression";

    // --- read counter qubits ---
    read_reg_int(istr, ids_counter);

    // --- Read the name of the circuit U whose product the gadget will compute ---
    istr >> name_U;
    if(ocs.find(name_U) == ocs.end())
        throw string("CompressionGadget: a circuit with the name ["s + name_U + "] is not found."s);
    YSQ oc_U = ocs[name_U];

    // --- Read qubits where the circuit U sits ---
    read_reg_int(istr, ids_U_target);

    // --- Read the integer indicating how many copies of U are multiplied ---
    try
    {
        istr >> word;
        N_mult = get_value_from_word(word);
    }
    catch(YCS e)
    {
        throw "CompressionGadget: wrong format of the N_mult.";
    }
    data.N_mult = N_mult;

    // --- Check whether the counter register has enough qubits ---
    int nc = ids_counter.size();
    if(N_mult == 1)
        if(nc < 1)
            throw string("CompressionGadget: if N_mult >= 1, "s + 
                "then n of qubits in the counter register has to be >= 1."s);
    if(N_mult == 2)
        if(nc < 2)
            throw string("CompressionGadget: if N_mult >= 2, "s + 
                "then n of qubits in the counter register has to be >= 2."s);
    // if(N_mult >= 3)
    // {
    //     int temp = ceil(log2(N_mult));
    //     temp = temp + (1 - ceil( (N_mult%int(pow(2,temp))) / N_mult)) + 1;
    //     if(nc < temp)
    //     {
    //         throw string("CompressionGadget: the number of qubits" + 
    //             " in the counter register is not large enough: "s + 
    //             "(nc-minimum = " + to_string(temp) + ").");
    //     }
    // }
        
    // --- read the flag whether output state after each call to oc_U should be written ---
    istr >> word;
    flag_step_output = get_value_from_word(word);

    // --- read the end of the gate structure ---
    read_end_gate(istr, ids_unit, ids_zero);

    // --- Print data ---
    cout << "\n--- Parameters of the compression gadget " << gadget_name << " ---" << endl;
    cout << "N_mult: " << data.N_mult << endl;
    cout << "Operator: " << name_U << endl;
    cout << "flag_Step_output: " << flag_step_output << endl;
    cout << "\n";

    // --- Construct the compression gadget ---
    compression_gadget(
        data, ids_counter, oc_U, ids_U_target, N_mult, flag_step_output, 
        ids_unit, ids_zero, flag_inv
    );

    // store the gadget data:
    if(!flag_data_read)
        map_gadget_data[gadget_name] = data;
}




void QCircuit::read_structure_repeat(YISS istr, std::map<std::string, YSQ>& ocs, YCB flag_inv)
{
    YVIv ids_U_target;
    string name_U;
    int N_mult;
    YVIv ids_unit, ids_zero; 
    string word;

    // --- Read the name of the circuit U to repeat ---
    istr >> name_U;
    if(ocs.find(name_U) == ocs.end())
        throw string("CompressionGadget: a circuit with the name ["s + name_U + "] is not found."s);
    YSQ oc_U = ocs[name_U];

    // --- Read qubits where the circuit U sits ---
    read_reg_int(istr, ids_U_target);

    // --- Read the integer indicating how many copies of U are multiplied ---
    try
    {
        istr >> word;
        N_mult = get_value_from_word(word);
    }
    catch(YCS e)
    {
        throw "CompressionGadget: wrong format of the N_mult.";
    }

    // --- read the end of the gate structure ---
    read_end_gate(istr, ids_unit, ids_zero);

    // --- Construct the compression gadget ---
    repeat(oc_U, ids_U_target, N_mult, ids_unit, ids_zero, flag_inv);
}






void QCircuit::read_structure_gate_phase_estimation(YISS istr, YCS path_in, std::map<std::string, YSQ>& ocs, YCB flag_inv)
{
    YVIv ids_ta; // target qubits of the operator A, whose eigenphase we seek for;
    YVIv ids_ty;
    YVIv ids_unit, ids_zero; 
    string name_A, name_INIT;

    // --- read target qubits of A ---
    read_reg_int(istr, ids_ta);

    // --- Read names of the operator A and INIT to find them among available circuits "ocs" ---
    try
    {
        istr >> name_A >> name_INIT;
    }
    catch(YCS e)
    {
        throw "PE definition: wrong format of the circuit names of the operators A and INIT: " + e;
    }

    // --- Read target qubit where the phase is to be written to ---
    read_reg_int(istr, ids_ty);
    
    // --- read the end of the gate structure ---
    read_end_gate(istr, ids_unit, ids_zero);

    // --- Find the appropriate circuits ---
    if(ocs.find(name_A) == ocs.end())
        throw string("PE definition: a circuit with the name ["s + name_A + "] is not found."s);
    if(ocs.find(name_INIT) == ocs.end())
        throw string("PE definition: a circuit with the name ["s + name_INIT + "] is not found."s);
    YSQ oc_A = ocs[name_A];
    YSQ oc_INIT = ocs[name_INIT];

    // --- add the phase estimation circuit ---
    phase_estimation(ids_ta, oc_A, oc_INIT, ids_ty, ids_unit, ids_zero, flag_inv);
}


void QCircuit::read_structure_gate_qsvt(
    YISS istr, 
    std::map<std::string, YSQ>& ocs, 
    YCB flag_inv, 
    QuCF_complex_data& qucf_d,
    YCB flag_QETU
){
    string name_circuit;
    string be_name;
    YVIv ids_a_qsvt, ids_be;
    YVIv ids_unit, ids_zero;
    auto& map_qsvt_data = qucf_d.qsvt;
    QSVT_pars data;

    bool flag_data_read;

    // --- read QSVT parameters ---
    istr >> name_circuit;
    flag_data_read = qucf_d.check_name(name_circuit);
    data.name = name_circuit;

    qsvt_read_parameters(name_circuit, data);

    // --- read QSVT ancilla ---
    read_reg_int(istr, ids_a_qsvt);
    if(data.parity >= 0)
    {
        if(ids_a_qsvt.size() != 1)
            throw "QSVT circuit of type ["s + data.type + "] must have only a single ancilla specific qubit."s;
    }
    else
    {
        if(ids_a_qsvt.size() != 2)
            throw "QSVT/QSP circuit of type ["s + data.type + "] must have two ancilla specific qubits."s;
    }
    
    // --- read the name of the block-encoding oracle ---
    istr >> be_name;
    if(ocs.find(be_name) == ocs.end())
        throw string("QSVT definition: a circuit with the name ["s + be_name + "] is not found."s);
    YSQ oc_be = make_shared<QCircuit>(ocs[be_name]);

    // auto oc_be_inv = make_shared<QCircuit>(oc_be);
    // oc_be_inv->h_adjoint();

    // --- read qubits where the BE oracle must be placed ---
    read_reg_int(istr, ids_be);
    if(oc_be->get_n_qubits() != ids_be.size())
    {
        stringstream temp_sstr;
        temp_sstr << "QSVT definition: Number of qubits in the oracle is " << oc_be->get_n_qubits() << ", ";
        temp_sstr << "but it is indicated that the oracle is placed on " << ids_be.size() << " qubits.\n";
        throw string(temp_sstr.str());
    }

    // --- read the end of the gate structure ---
    read_end_gate(istr, ids_unit, ids_zero);

    // --- add the QSVT circuit ---
    if(flag_QETU)
    {
        cout << "Constructing QETU circuit..." << endl;
        vector<double> angles_qetu;
        if(data.parity == 0)
            angles_qetu = data.angles_phis_even;
        if(data.parity == 1)
            angles_qetu = data.angles_phis_odd;

        // // correct angles: from QSVT to QETU:
        // int Na = angles_qetu.size();
        // angles_qetu[0]    += M_PI_4;
        // angles_qetu[Na-1] += M_PI_4;
        // for(int ii = 1; ii < (Na-1); ii++)
        //     angles_qetu[ii] += M_PI_2;
 
        // construct the QET circuit: 
        qetu_def_parity(angles_qetu, ids_a_qsvt[0], ids_be, oc_be, ids_unit, ids_zero, flag_inv);
    }
    else
    {
        cout << "Constructing QSVT circuit..." << endl;
        if(data.parity == 0)
            qsvt_def_parity(data.angles_phis_even, ids_a_qsvt[0], ids_be, oc_be, ids_unit, ids_zero, flag_inv);
        if(data.parity == 1)
            qsvt_def_parity(data.angles_phis_odd, ids_a_qsvt[0], ids_be, oc_be, ids_unit, ids_zero, flag_inv); 
        if(data.parity == -1)
        {
            // cout << "a-qsp: " << ids_a_qsvt[1] << endl;
            // cout << "a-qu: "  << ids_a_qsvt[0] << endl;
            // cout << "be[0]: " << ids_be[0] << endl;
            if(YMIX::compare_strings(data.type, "QSP-ham")){
                // qsp_ham(
                //     name_circuit,
                //     data.angles_phis_arbitrary, 
                //     data.n_repeat,
                //     ids_a_qsvt[1], ids_a_qsvt[0], 
                //     ids_be, 
                //     oc_be, 
                //     ids_unit, ids_zero, 
                //     flag_inv
                // );
                qsp_ham_opt2(
                    name_circuit,
                    data.angles_phis_arbitrary, 
                    data.n_repeat,
                    ids_a_qsvt[1], ids_a_qsvt[0], 
                    ids_be, 
                    oc_be, 
                    ids_unit, ids_zero, 
                    flag_inv
                );
            }
        }
    }
    

    // store the QSVT data:
    if(!flag_data_read)
        map_qsvt_data[name_circuit] = data;
}


void QCircuit::read_selector_power(YISS istr, std::map<std::string, YSQ>& ocs, YCB flag_inv)
{
    YVIv rs, ids_U_target, ids_unit, ids_zero;
    string name_U;

    // --- read selector qubits ---
    read_reg_int(istr, rs);

    // --- Read the name of the circuit U whose powers will be computed ---
    istr >> name_U;
    if(ocs.find(name_U) == ocs.end())
        throw string("SelectorPower: a circuit with the name ["s + name_U + "] is not found."s);
    YSQ oc_U = ocs[name_U];

    // --- Read qubits where the circuit U sits ---
    read_reg_int(istr, ids_U_target);

    // --- read the end of the gate structure ---
    read_end_gate(istr, ids_unit, ids_zero);

    // --- Construct the compression gadget ---
    selector_power(rs, oc_U, ids_U_target, ids_unit, ids_zero, flag_inv);
}


void QCircuit::read_structure_LCHS_QSP(YISS istr, std::map<std::string, YSQ>& ocs, YCB flag_inv)
{
    YVIv rs, ids_Ut, ids_Ow, ids_unit, ids_zero;
    string name_U, word;
    int Nt;
    YSQ Ut, Ow;

    // --- Read the name of the oracle simulating a short time interval ---
    istr >> name_U;
    if(ocs.find(name_U) == ocs.end())
        throw string("LCHS-QSP: a circuit with the name ["s + name_U + "] is not found."s);
    Ut = ocs[name_U];

    // --- Read qubits where the circuit Ut sits ---
    read_reg_int(istr, ids_Ut);

    // --- Read the number of time steps ---
    istr >> word;
    Nt = get_value_from_word(word);

    // --- Read the name of the oracle computing weights ---
    istr >> name_U;
    if(ocs.find(name_U) == ocs.end())
        throw string("LCHS-QSP: a circuit with the name ["s + name_U + "] is not found."s);
    Ow = ocs[name_U]; 

    // --- Read qubits where the circuit Ow sits ---
    read_reg_int(istr, ids_Ow);

    // --- read the end of the gate structure ---
    read_end_gate(istr, ids_unit, ids_zero);

    // --- Construct the LCHS-QSP circuit ---
    LCHS_QSP(
        Ut, ids_Ut, Nt,
        Ow, ids_Ow,
        ids_unit, ids_zero, flag_inv
    );
}


void QCircuit::read_structure_dirdec(YISS istr, YCB flag_inv)
{
    using namespace std::complex_literals;

    YVIv ids_ctrl, cs_unit, cs_zero;
    int id_targ;
    string sel_prof, word;
    int N_pars;
    YVQv pars;

    // --- Read the target qubit ---
    YVIv temp_vec;
    read_reg_int(istr, temp_vec);
    id_targ = temp_vec[0];

    // --- Read control qubits ---
    read_reg_int(istr, ids_ctrl);

    // --- Read the string defining the profile to compute ---
    istr >> sel_prof;
    YMIX::print_log("DirDec: compute the profile: " + sel_prof);

    // --- Read the number of parameters for the target profile ---
    istr >> word;
    N_pars = get_value_from_word(word);
    YMIX::print_log("DirDec: N_pars = " + to_string(N_pars));

    // --- Read the parameters for the target profile ---
    double v1;
    for(int i_par = 0; i_par < N_pars; i_par++)
    {
        istr >> word;
        v1 = get_value_from_word(word);
        pars.push_back(v1);
        YMIX::print_log("DirDec: par" + to_string(i_par) + " = " + to_string(v1));
    }

    // --- read the end of the gate structure ---
    read_end_gate(istr, cs_unit, cs_zero);

    // --- Construct direct decomposition ---
    int n_ctrl = ids_ctrl.size();
    int Nc = 1 << n_ctrl;
    if(YMIX::compare_strings(sel_prof, "LCHS_weights_sqrt"))
    {
        if(N_pars < 1)
            throw std::string("ERROR: DirDec LCHS_weights_sqrt: N_pars should be = 1");

        YVQv phis_y(Nc);
        double k;
        double dk = 2./(Nc - 1);
        for(int ii = 0; ii < Nc; ii++)
        {
            k = -1.0 + ii * dk;
            k *= pars[0];
            v1 = 1. / sqrt(1 + k*k);
            phis_y[ii] = 2. * acos(v1);
        }
        DirDec_Y(phis_y, ids_ctrl, id_targ, cs_unit, cs_zero, flag_inv);
    }
    else if(YMIX::compare_strings(sel_prof, "LCHS_weights_sin_sqrt"))
    {
        if(N_pars < 1)
            throw std::string("ERROR: DirDec LCHS_weights_sin_sqrt: N_pars should be = 1");

        YVQv phis_y(Nc);
        double theta_1, k;
        double d_theta = M_PI/(Nc - 1);
        double kmax = pars[0];
        for(int ii = 0; ii < Nc; ii++)
        {
            theta_1 = -M_PI/2. + ii * d_theta;
            k  = kmax * sin(theta_1);
            v1 = kmax * cos(theta_1) / (1 + k*k);
            v1 = sqrt(v1 * d_theta / M_PI);
            phis_y[ii] = 2. * acos(v1);
        }
        DirDec_Y(phis_y, ids_ctrl, id_targ, cs_unit, cs_zero, flag_inv);
    }
    else if(YMIX::compare_strings(sel_prof, "LCHS_weights_sin_OPT_sqrt"))
    {
        if(N_pars < 2)
            throw std::string("ERROR: DirDec LCHS_weights_sin_OPT_sqrt: N_pars should be = 2");

        YVQv phis_y(Nc);
        YVQv phis_z(Nc);
        double theta_1, k;
        double d_theta = M_PI/(Nc - 1);
        double kmax = pars[0];
        double beta = pars[1];
        complex<double> c1, c2;
        double coef_beta = 2. * M_PI * exp(-pow(2,beta));
        for(int ii = 0; ii < Nc; ii++)
        {
            theta_1 = -M_PI/2. + ii * d_theta;
            k  = kmax * sin(theta_1);
            c1 = kmax * cos(theta_1) * d_theta;
            c2 = coef_beta * exp(pow(1.+1i*k, beta)) * (1. - 1i * k); 
            c1 = sqrt(c1 / c2);
            // c1 = c1 / c2;
            // cout << "ii, w: " << ii << c1 << endl;
            phis_y[ii] = 2. * acos(abs(c1));
            phis_z[ii] = -2. * arg(c1);
        }
        DirDec_C(phis_y, phis_z, ids_ctrl, id_targ, cs_unit, cs_zero, flag_inv);
    }
    else if(YMIX::compare_strings(sel_prof, "LCHS_weights_full"))
    {
        if(N_pars < 1)
            throw std::string("ERROR: DirDec LCHS_weights_full: N_pars should be = 1");

        YVQv phis_y(Nc);
        double k;
        double dk = 2./(Nc - 1);
        for(int ii = 0; ii < Nc; ii++)
        {
            k = -1.0 + ii * dk;
            k *= pars[0];
            v1 = 1. / (1 + k*k);
            phis_y[ii] = 2. * acos(v1);
        }
        DirDec_Y(phis_y, ids_ctrl, id_targ, cs_unit, cs_zero, flag_inv);
    }
    else if(YMIX::compare_strings(sel_prof, "linear"))
    {
        if(N_pars < 2)
            throw std::string("ERROR: DirDec linear: N_pars should be = 2");

        YVQv phis_y(Nc);
        double k;
        double dk = (pars[1] - pars[0])/(Nc - 1);
        for(int ii = 0; ii < Nc; ii++)
        {
            k = pars[0] + ii * dk;
            phis_y[ii] = 2. * acos(k);
        }
        DirDec_Y(phis_y, ids_ctrl, id_targ, cs_unit, cs_zero, flag_inv);
    }
    else if(YMIX::compare_strings(sel_prof, "gauss"))
    {
        // gauss in the interval [0, 1];
        if(N_pars < 2)
            throw std::string("ERROR: DirDec gauss: N_pars should be = 2");

        YVQv phis_y(Nc);
        double xc = pars[0];
        double wd = pars[1];
        double x, f;
        double dx = 1./(Nc - 1);
        for(int ii = 0; ii < Nc; ii++)
        {
            x = 0 + ii * dx;
            f = (x - xc) * (x - xc) / (2. * wd*wd);
            f = exp(-f);
            phis_y[ii] = 2. * acos(f);
        }
        DirDec_Y(phis_y, ids_ctrl, id_targ, cs_unit, cs_zero, flag_inv);
    }
    else if(YMIX::compare_strings(sel_prof, "inverse"))
    {
        if(N_pars < 3)
            throw std::string("ERROR: DirDec inverse: N_pars should be = 3");

        YVQv phis_y(Nc);
        int Nch = 1 << (n_ctrl - 1);
        double kappa = pars[0];
        double xmax  = pars[1];
        double norm_coef = pars[2];
        double x0 = 1./kappa;
        double dx = (xmax - x0) / (Nch - 1);
        double x1;
        for(int ii = 0; ii < Nc; ii++)
        {
            if(ii < Nch)
                x1 = -xmax + dx * ii;
            else
                x1 = x0 + dx * (ii - Nch);
            v1 = 1. - exp(-pow(5.*kappa*x1, 2.));
            v1 /= (x1 * kappa);
            v1 *= norm_coef;
            phis_y[ii] = 2. * acos(v1);
        }
        DirDec_Y(phis_y, ids_ctrl, id_targ, cs_unit, cs_zero, flag_inv);
    }
    else
    {
        YMIX::print_log(
            ">>> WARNING: DirDec: the profile " + sel_prof + 
            " is not recognized. Skipping the DirDec."
        );  
    }
}




YQCP QCircuit::adder_1(YCVI ts, YCVI cs_unit, YCVI cs_zero, YCB flag_inv)
{
    YVIv ids_target = YVIv(ts);
    uint32_t nt;

    if(!flag_inv)
    {
        // put the high-priority qubits at the beginning
        sort(ids_target.begin(), ids_target.end(), greater<int>());
        nt = ids_target.size();

        // add CNOT and X gates with control nodes
        for(unsigned i = 0; i < nt-1; ++i)
        {
            YVIv ids_cnot_cs = YVIv(ids_target.begin() + i + 1, ids_target.end());
            ids_cnot_cs.insert(ids_cnot_cs.end(), cs_unit.begin(), cs_unit.end());
            x(ids_target[i], ids_cnot_cs, cs_zero);
        }
        x(ids_target.back(), cs_unit, cs_zero);
    }
    else
    {
        subtractor_1(ts, cs_unit, cs_zero, false);
    }
    return get_the_circuit();
}


YQCP QCircuit::adder_2(YCVI ts, YCVI cs_unit, YCVI cs_zero, YCB flag_inv)
{
    YVIv ids_target_init = YVIv(ts);
    uint32_t nt;

    if(!flag_inv)
    {
        // put the high-priority qubits at the beginning
        sort(ids_target_init.begin(), ids_target_init.end(), greater<int>());
        nt = ids_target_init.size()-1;
        YVIv ids_target(nt);
        copy(ids_target_init.begin(), ids_target_init.end()-1, ids_target.begin());
        
        // add CNOT and X gates with control nodes
        for(unsigned i = 0; i < nt-1; ++i)
        {
            YVIv ids_cnot_cs = YVIv(ids_target.begin() + i + 1, ids_target.end());
            ids_cnot_cs.insert(ids_cnot_cs.end(), cs_unit.begin(), cs_unit.end());
            x(ids_target[i], ids_cnot_cs, cs_zero);
        }
        x(ids_target.back(), cs_unit, cs_zero);
    }
    else
    {
        subtractor_2(ts, cs_unit, cs_zero, false);
    }
    return get_the_circuit();
}


YQCP QCircuit::adder_3(YCVI ts, YCVI cs_unit, YCVI cs_zero, YCB flag_inv)
{
    YVIv ids_target = YVIv(ts);
    int32_t nt;
    int32_t sh = 2;

    if(!flag_inv)
    {
        // put the high-priority qubits at the beginning
        sort(ids_target.begin(), ids_target.end(), greater<int>());
        nt = ids_target.size();

        // add CNOT and X gates with control nodes
        for(int i = 0; i < (nt-1-sh); ++i)
        {
            YVIv ids_cnot_cs = YVIv(ids_target.begin() + i + 1, ids_target.end()-sh);
            ids_cnot_cs.insert(ids_cnot_cs.end(), cs_unit.begin(), cs_unit.end());
            x(ids_target[i], ids_cnot_cs, cs_zero);
        }
        if((nt-1-sh) >= 0) 
            x(ids_target[nt-1-sh], cs_unit, cs_zero);

        x(ids_target.back(), cs_unit, cs_zero);
        for(int i = nt-2; i >= 0; --i)
        {
            YVIv ids_cnot_cs = YVIv(ids_target.begin() + i + 1, ids_target.end());
            ids_cnot_cs.insert(ids_cnot_cs.end(), cs_unit.begin(), cs_unit.end());
            x(ids_target[i], ids_cnot_cs, cs_zero);
        }
    }
    else
    {
        subtractor_3(ts, cs_unit, cs_zero, false);
    }
    return get_the_circuit();
}


YQCP QCircuit::subtractor_1(YCVI ts, YCVI cs_unit, YCVI cs_zero, YCB flag_inv)
{
    YVIv ids_target = YVIv(ts);
    uint32_t nt;

    if(!flag_inv)
    {
        // put the low-priority qubits at the beginning 
        sort(ids_target.begin(), ids_target.end());
        nt = ids_target.size();

        // add CNOT and X gates with control nodes 
        x(ids_target[0], cs_unit, cs_zero);
        for(unsigned i = 1; i < nt; ++i)
        {
            YVIv ids_cnot_cs = YVIv(ids_target.begin(), ids_target.begin() + i);
            ids_cnot_cs.insert(ids_cnot_cs.end(), cs_unit.begin(), cs_unit.end());
            x(ids_target[i], ids_cnot_cs, cs_zero);
        }
    }
    else
    {
        adder_1(ts, cs_unit, cs_zero, false);
    }
    return get_the_circuit();
}


YQCP QCircuit::subtractor_2(YCVI ts, YCVI cs_unit, YCVI cs_zero, YCB flag_inv)
{
    YVIv ids_target_init = YVIv(ts);
    uint32_t nt;

    if(!flag_inv)
    {
        // put the low-priority qubits at the beginning 
        sort(ids_target_init.begin(), ids_target_init.end());
        nt = ids_target_init.size() - 1;
        YVIv ids_target(nt);
        copy(ids_target_init.begin()+1, ids_target_init.end(), ids_target.begin());

        // add CNOT and X gates with control nodes 
        x(ids_target[0], cs_unit, cs_zero);
        for(unsigned i = 1; i < nt; ++i)
        {
            YVIv ids_cnot_cs = YVIv(ids_target.begin(), ids_target.begin() + i);
            ids_cnot_cs.insert(ids_cnot_cs.end(), cs_unit.begin(), cs_unit.end());
            x(ids_target[i], ids_cnot_cs, cs_zero);
        }
    }
    else
    {
        adder_2(ts, cs_unit, cs_zero, false);
    }
    return get_the_circuit();
}


YQCP QCircuit::subtractor_3(YCVI ts, YCVI cs_unit, YCVI cs_zero, YCB flag_inv)
{
    YVIv ids_target = YVIv(ts);
    uint32_t nt;
    uint32_t sh = 2;

    if(!flag_inv)
    {
        // put the low-priority qubits at the beginning 
        sort(ids_target.begin(), ids_target.end());
        nt = ids_target.size();

        for(unsigned i = nt-1; i > 0; --i)
        {
            YVIv ids_cnot_cs = YVIv(ids_target.begin(), ids_target.begin() + i);
            ids_cnot_cs.insert(ids_cnot_cs.end(), cs_unit.begin(), cs_unit.end());
            x(ids_target[i], ids_cnot_cs, cs_zero);
        }
        x(ids_target[0], cs_unit, cs_zero);

        if(nt > sh) x(ids_target[sh], cs_unit, cs_zero);
        for(unsigned i = sh+1; i < nt; ++i)
        {
            YVIv ids_cnot_cs = YVIv(ids_target.begin() + sh, ids_target.begin() + i);
            ids_cnot_cs.insert(ids_cnot_cs.end(), cs_unit.begin(), cs_unit.end());
            x(ids_target[i], ids_cnot_cs, cs_zero);
        }
    }
    else
    {
        adder_3(ts, cs_unit, cs_zero, false);
    }
    return get_the_circuit();
}


YQCP QCircuit::adder(YCVI ts1, YCVI ts2, YCVI ts3, YCVI cs_unit, YCVI cs_zero, YCB flag_inv, YCB flag_box)
{
    uint32_t nt = ts1.size();
    string oracle_name_tex = "ADD";
    YVIv ts_total; 
    ts_total.insert(ts_total.end(), ts1.begin(), ts1.end());
    ts_total.insert(ts_total.end(), ts2.begin(), ts2.end());
    ts_total.insert(ts_total.end(), ts3.begin()+1, ts3.end());
    ts_total.insert(ts_total.end(), ts3.begin(),   ts3.begin()+1);
    uint32_t nt_total = ts_total.size();

    // --- create an envelop circuit for the adder operator ---
    auto oc_adder = make_shared<QCircuit>(
        "ADD", env_, path_to_output_, nt_total
    );
    auto rc = oc_adder->add_register("rc", nt);
    auto rb = oc_adder->add_register("rb", nt);
    auto ra = oc_adder->add_register("ra", nt);

    for(uint32_t ii = 0; ii < nt; ii++)
    {
        oc_adder->x(rc[ii], YVIv {rb[ii], ra[ii]});
        oc_adder->x(rb[ii], YVIv {ra[ii]});
        if(ii > 0)
            oc_adder->x(rc[ii], YVIv {rc[ii-1], rb[ii]});
    } 

    if(nt > 1) 
    {
        oc_adder->x(rb[nt-1], YVIv {rc[nt-2]});
        for(uint32_t ii = nt-2; ii > 0; ii--)
        {
            oc_adder->x(rc[ii], YVIv {rc[ii-1], rb[ii]});
            oc_adder->x(rb[ii], YVIv {ra[ii]});
            oc_adder->x(rc[ii], YVIv {rb[ii], ra[ii]});
            oc_adder->x(rb[ii], YVIv {ra[ii]});
            oc_adder->x(rb[ii], YVIv {rc[ii-1]});
        }
        oc_adder->x(rb[0], YVIv {ra[0]});
        oc_adder->x(rc[0], YVIv {rb[0], ra[0]});
        oc_adder->x(rb[0], YVIv {ra[0]});
    }

    // --- invert the circuit if necessary ---
    if(flag_inv)
    {
        oc_adder->h_adjoint();
        oracle_name_tex += "^\\dagger";
    }

    // --- copy the adder circuit to the current circuit ---
    auto box = YSB(nullptr);
    if(flag_box)
        box = YMBo("ADD", ts_total, YVIv{}, YVIv{}, oracle_name_tex);
    copy_gates_from(
        oc_adder,
        ts_total,
        box, 
        cs_unit, cs_zero,
        false         
    );
    return get_the_circuit();
}


YQCP QCircuit::subtractor(YCVI ts1, YCVI ts2, YCVI ts3, YCVI cs_unit, YCVI cs_zero, YCB flag_inv)
{
    x(ts2, cs_unit, cs_zero);
    adder(ts1, ts2, ts3, cs_unit, cs_zero, flag_inv);
    x(ts2, cs_unit, cs_zero);
    return get_the_circuit();
}


YQCP QCircuit::adder_qft(
    YCVI qs_v1, YCVI qs_v2, YCI id_carry, 
    YCVI cs_unit, YCVI cs_zero, 
    YCB flag_inv, YCB flag_box
){
    int nv       = qs_v1.size();
    int nt_total = nv + 1;
    int n_circ_total = nt_total + nv;
    string oracle_name_tex = "ADDQFT";

    // --- create an envelop circuit for the adder operator ---
    auto oc_adder = make_shared<QCircuit>(
        "ADD-QFT", env_, path_to_output_, n_circ_total
    );
    auto loc_v2    = oc_adder->add_register("a", nv);
    auto loc_carry = oc_adder->add_register("c", 1)[0];
    auto loc_v1    = oc_adder->add_register("r", nv);
    auto loc_res = YVIv(loc_v1);
    loc_res.push_back(loc_carry);

    oc_adder->quantum_fourier(loc_res, YVIv{}, YVIv{}, false, true);

    // starting from the most significant qubit (which is the carry bit):
    int n_phase_gates_per_qubit = 0;
    for(int i_qubit = nt_total-1; i_qubit >= 0; i_qubit--)
    {
        n_phase_gates_per_qubit++;
        for(int j_ph_gate = 0; j_ph_gate < n_phase_gates_per_qubit; j_ph_gate++)
        {
            if(j_ph_gate < nv)
            {
                double phase_curr = 2*M_PI / (1<<(n_phase_gates_per_qubit - j_ph_gate));
                oc_adder->phase(i_qubit, phase_curr, YVIv{loc_v2[j_ph_gate]});
            }
        }
    }
    oc_adder->quantum_fourier(loc_res, YVIv{}, YVIv{}, true, true);

    // --- invert the circuit if necessary ---
    if(flag_inv)
    {
        oc_adder->h_adjoint();
        oracle_name_tex += "^\\dagger";
    }

    // --- copy the adder circuit to the current circuit ---
    auto ts_total = YVIv(qs_v1);
    ts_total.push_back(id_carry);
    ts_total.insert(ts_total.end(), qs_v2.begin(), qs_v2.end());

    auto box = YSB(nullptr);
    if(flag_box)
        box = YMBo("ADDFIXED", ts_total, YVIv{}, YVIv{}, oracle_name_tex);
    copy_gates_from(
        oc_adder,
        ts_total,
        box, 
        cs_unit, cs_zero,
        false   
    );
    return get_the_circuit();
}


YQCP QCircuit::subtractor_qft(
    YCVI qs_v1, YCVI qs_v2, YCI id_sign, 
    YCVI cs_unit, YCVI cs_zero, 
    YCB flag_inv, YCB flag_box
){

    adder_qft(qs_v1, qs_v2, id_sign, cs_unit, cs_zero, !flag_inv, flag_box);

    // x(qs_v1, cs_unit, cs_zero);
    // adder_qft(qs_v1, qs_v2, id_sign, cs_unit, cs_zero, flag_inv, flag_box);
    // x(qs_v1, cs_unit, cs_zero);
    return get_the_circuit();
}


YQCP QCircuit::adder_fixed(
    YCVI ids_target, YCI id_carry, YCU int_sub, 
    YCVI cs_unit, YCVI cs_zero, 
    YCB flag_inv, YCB flag_box
){
    int nv       = ids_target.size();
    int nt_total = nv + 1;
    string oracle_name_tex = "ADDFIXED";

    // --- represent the unisgned integer as a bitstring ---
    YVshv bitstring_int(nt_total);
    YMATH::intToBinary(int_sub, bitstring_int);

    // --- create an envelop circuit for the adder operator ---
    auto oc_adder = make_shared<QCircuit>(
        "ADD-FIXED", env_, path_to_output_, nt_total
    );
    auto r_all = oc_adder->add_register("r", nt_total);

    oc_adder->quantum_fourier(r_all);

    // starting from the most significant qubit (which is the carry bit):
    int n_phase_gates_per_qubit = 0;
    for(int i_qubit = nt_total-1; i_qubit >= 0; i_qubit--)
    {
        n_phase_gates_per_qubit++;
        for(int j_ph_gate = 0; j_ph_gate < n_phase_gates_per_qubit; j_ph_gate++)
            if(bitstring_int[nt_total - j_ph_gate - 1] == 1)
            {
                double phase_curr = 2*M_PI / (1<<(n_phase_gates_per_qubit - j_ph_gate));
                oc_adder->phase(i_qubit, phase_curr);
            }
    }
    oc_adder->quantum_fourier(r_all, YVIv{}, YVIv{}, true);

    // --- invert the circuit if necessary ---
    if(flag_inv)
    {
        oc_adder->h_adjoint();
        oracle_name_tex += "^\\dagger";
    }

    // --- copy the adder circuit to the current circuit ---
    auto ts_total = YVIv(ids_target);
    ts_total.push_back(id_carry);

    auto box = YSB(nullptr);
    if(flag_box)
        box = YMBo("ADDFIXED", ts_total, YVIv{}, YVIv{}, oracle_name_tex);
    copy_gates_from(
        oc_adder,
        ts_total,
        box, 
        cs_unit, cs_zero,
        false   
    );
    return get_the_circuit();
}


YQCP QCircuit::subtractor_fixed(
        YCVI ids_target, YCI id_carry, YCU int_sub, 
        YCVI cs_unit, YCVI cs_zero,
        YCB flag_inv
){
    auto ts_total = YVIv(ids_target);
    ts_total.push_back(id_carry);        
    // x(ts_total, cs_unit, cs_zero);
    // adder_fixed(ids_target, id_carry, int_sub, cs_unit, cs_zero, flag_inv);
    // x(ts_total, cs_unit, cs_zero);

    adder_fixed(ids_target, id_carry, int_sub, cs_unit, cs_zero, !flag_inv);
    return get_the_circuit();
}


YQCP QCircuit::comparator_fixed(
        YCVI ids_target, YCVI ids_carry, YCU int_sub, 
        YCVI cs_unit, YCVI cs_zero,
        YCB flag_inv
){
    auto cs_total_for_x = YVIv(cs_unit);
    cs_total_for_x.push_back(ids_carry[0]);

    subtractor_fixed(ids_target, ids_carry[0], int_sub, cs_unit, cs_zero);
    x(ids_carry[1], cs_total_for_x, cs_zero);
    subtractor_fixed(ids_target, ids_carry[0], int_sub, cs_unit, cs_zero, true);
    return get_the_circuit();
}


YQCP QCircuit::swap(YCI t1, YCI t2, YCVI cs_unit, YCVI cs_zero)
{ 
    YVIv ids_cs_1_unit = {t1};
    YVIv ids_cs_2_unit = {t2};
    ids_cs_1_unit.insert(ids_cs_1_unit.end(), cs_unit.begin(), cs_unit.end());
    ids_cs_2_unit.insert(ids_cs_2_unit.end(), cs_unit.begin(), cs_unit.end());
    return x(t2, ids_cs_1_unit, cs_zero)->x(t1, ids_cs_2_unit, cs_zero)->x(t2, ids_cs_1_unit, cs_zero); 
}


YQCP QCircuit::quantum_fourier(YCVI ts, YCVI cs_unit, YCVI cs_zero, YCB flag_inv, YCB flag_box)
{
    uint32_t nt = ts.size();
    uint32_t q1;
    int cq;
    qreal aa;
    string fourier_name_tex = "F";

    // --- create an envelop circuit for the Fourier operator ---
    auto oc_fourier = make_shared<QCircuit>(
        "F", env_, path_to_output_, nt
    );
    auto q = oc_fourier->add_register("q", nt);
    for(uint32_t ii = 0; ii < nt; ii++)
    {
        q1 = nt - 1 - ii;
        oc_fourier->h(q1);
        for(uint32_t jj = 1; jj < (nt - ii); jj++)
        {
            aa = M_PI / (1 << jj);
            oc_fourier->phase(q1, aa, YVIv{int(q1 - jj)});
        }
    }
    for(uint32_t ii = 0; ii < uint32_t(nt/2); ii++)
        oc_fourier->swap(q[ii], q[nt - 1 - ii]);

    // --- invert the circuit if necessary ---
    if(flag_inv)
    {
        oc_fourier->h_adjoint();
        fourier_name_tex += "^\\dagger";
    }

    // --- copy Fourier env. circuit to the current circuit ---
    auto box = YSB(nullptr);
    if(flag_box)
        box = YMBo("F", ts, YVIv{}, YVIv{}, fourier_name_tex);
    copy_gates_from(
        oc_fourier,
        ts,
        box, 
        cs_unit, cs_zero,
        false         
    );
    return get_the_circuit();
}


YQCP QCircuit::gate_sin(
        YCVI anc, 
        YCVI conds, 
        YCQR alpha_0, 
        YCQR alpha, 
        YCVI cs_unit, YCVI cs_zero, 
        YCB flag_inv, 
        YCB flag_box
){
    auto na = 1;
    auto n_cond = conds.size();
    auto n_tot = n_cond + na;
    string name_tex = "SIN";

    // --- all target qubits ---
    vector<int> qubits_tot = YVIv(conds);
    qubits_tot.insert(qubits_tot.end(), anc.begin(), anc.end());

    // --- create an envelop circuit for the sin gate ---
    auto oc_sin = make_shared<QCircuit>(
        "SIN", env_, path_to_output_, n_tot
    );
    auto a_loc    = oc_sin->add_register("a", na)[0];
    auto cond_loc = oc_sin->add_register("cond", n_cond);

    oc_sin->ry(a_loc, 2*alpha_0);
    for(auto ii = 0; ii < n_cond; ii++)
    {
        qreal aa = 2*alpha / pow(2., n_cond - 1 - ii);
        oc_sin->ry(a_loc, aa, YVIv{cond_loc[ii]});
    }
    oc_sin->x(a_loc);

    // --- invert the circuit if necessary ---
    if(flag_inv)
    {
        oc_sin->h_adjoint();
        name_tex += "^\\dagger";
    }

    // --- copy the env. circuit to the current circuit ---
    auto box = YSB(nullptr);
    if(flag_box)
        box = YMBo("SIN", qubits_tot, YVIv{}, YVIv{}, name_tex);
    copy_gates_from(
        oc_sin,
        qubits_tot,
        box, 
        cs_unit, cs_zero,
        false        
    );
    return get_the_circuit();
}


YQCP QCircuit::gate_sinC(
        YCVI anc, 
        YCVI conds, 
        YCQR alpha_0_y, YCQR alpha_y, 
        YCQR alpha_0_z, YCQR alpha_z,
        YCVI cs_unit, YCVI cs_zero, YCB flag_inv
){
    auto na = 1;
    auto n_cond = conds.size();
    auto n_tot = n_cond + na;
    string name_tex = "SIN";

    // --- all target qubits ---
    vector<int> qubits_tot = YVIv(conds);
    qubits_tot.insert(qubits_tot.end(), anc.begin(), anc.end());

    // --- create an envelop circuit for the sin gate ---
    auto oc_sin = make_shared<QCircuit>(
        "SIN_C", env_, path_to_output_, n_tot
    );
    auto a_loc    = oc_sin->add_register("a", na)[0];
    auto cond_loc = oc_sin->add_register("cond", n_cond);

    // --- ampl: Ry rotations ---
    oc_sin->ry(a_loc, 2*alpha_0_y);
    for(auto ii = 0; ii < n_cond; ii++)
    {
        qreal aa = 2*alpha_y / pow(2., n_cond - 1 - ii);
        oc_sin->ry(a_loc, aa, YVIv{cond_loc[ii]});
    }

    // --- phase: Rz rotations ---
    oc_sin->rz(a_loc, 2*alpha_0_z);
    for(auto ii = 0; ii < n_cond; ii++)
    {
        qreal aa = 2*alpha_z / pow(2., n_cond - 1 - ii);
        oc_sin->rz(a_loc, aa, YVIv{cond_loc[ii]});
    }

    // --- choose sin(y) * exp(iz) ---
    oc_sin->x(a_loc);

    // --- invert the circuit if necessary ---
    if(flag_inv)
    {
        oc_sin->h_adjoint();
        name_tex += "^\\dagger";
    }

    // --- copy the env. circuit to the current circuit ---
    auto box = YSB(nullptr);
    copy_gates_from(
        oc_sin,
        qubits_tot,
        box, 
        cs_unit, cs_zero,
        false        
    );
    return get_the_circuit();
}




YQCP QCircuit::compression_gadget(
    GADGET_pars& data ,
    YCVI ids_counter, 
    YCCQ& oc_U_in, 
    YCVI ids_U_target, 
    YCI N_mult, 
    YCB flag_step_output,
    YCVI cs_unit, YCVI cs_zero,
    YCB flag_inv
){
    string gadget_name = data.name;
    string name_stop;
    auto oc_U = make_shared<QCircuit>(oc_U_in);
    auto oc_gadget = make_shared<QCircuit>("CG", env_, path_to_output_, nq_);
    oc_gadget->add_register("r", nq_);

    auto nq_U = oc_U->get_n_qubits();
    if(nq_U != ids_U_target.size())
    {
        throw string("Error in compression_gadget: the indicated number of target qubits for the oracle\n") +
            string("does not equal the number of qubits that the oracle requires.");
    }
    auto nanc_U = oc_U->get_na();

    YVIv inq_U(ids_U_target.begin(), ids_U_target.begin() + nq_U - nanc_U);
    YVIv anc_U(ids_U_target.begin() + nq_U - nanc_U, ids_U_target.end());

    // --- Apply the adder ---
    YVIv ids_t_adder;
    int id_t_carry;
    if(N_mult == 1 || N_mult == 2)
    {
        ids_t_adder = ids_counter;
        for(short i_inc = 0; i_inc < N_mult; i_inc++)
            oc_gadget->adder_1(ids_t_adder);
    }
    else
    {
        ids_t_adder = YVIv(ids_counter.begin(),   ids_counter.end()-1);
        id_t_carry  = ids_counter.back();
        oc_gadget->adder_fixed(ids_t_adder, id_t_carry, N_mult);
    }
    data.counter_qubits = ids_t_adder;
    
    // --- Calls to the operator oc_U and the subtractors ---
    for(int i_mult = 0; i_mult < N_mult; i_mult++)
    {
        oc_U->activate_gadget(N_mult-i_mult, N_mult);
        oc_gadget->copy_gates_from(oc_U, ids_U_target);
        oc_gadget->subtractor_1(ids_t_adder, YVIv{}, anc_U);
        if(flag_step_output && flag_stop_gates_)
        {
            name_stop = "CompressionGadget <" + gadget_name + ">: i = " + to_string(i_mult);
            oc_gadget->add_stop_gate(name_stop);
            // cout << "CG: add stop gate" << endl;
        }
    }

    // --- Transfer the gates to the main circuit ---
    auto all_qubits = YMATH::get_range(0, nq_);
    copy_gates_from(
        oc_gadget,
        all_qubits,
        YSB(nullptr), 
        cs_unit, cs_zero,
        flag_inv        
    ); 
    return get_the_circuit();
}



YQCP QCircuit::repeat(
    YCCQ& oc_U_in, 
    YCVI ids_U_target, 
    YCI N_mult, 
    YCVI cs_unit, YCVI cs_zero,
    YCB flag_inv
){
    YMIX::YTimer timer;
    timer.StartPrint("Creating the gate Repeat... ");
    if(flag_repeat_insert_)
    {
        if(flag_inv)
        {
            auto oc_U = make_shared<QCircuit>(oc_U_in);
            oc_U->h_adjoint();
            for(int i_mult = 0; i_mult < N_mult; i_mult++)
                insert_gates_from(oc_U.get());
        }
        else
        {
            for(int i_mult = 0; i_mult < N_mult; i_mult++)
                insert_gates_from(oc_U_in.get());
        }
    }
    else
    {
        auto nq_U = oc_U_in->get_n_qubits();
        if(nq_U != ids_U_target.size())
        {
            throw string("Error in repeat: the indicated number of target qubits for the oracle\n") +
                string("does not equal the number of qubits that the oracle requires.");
        }

        // auto oc_U = make_shared<QCircuit>(oc_U_in);
        // auto oc_repeat = make_shared<QCircuit>("CG", env_, path_to_output_, nq_);
        // oc_repeat->add_register("r", nq_);
        // for(int i_mult = 0; i_mult < N_mult; i_mult++)
        //     oc_repeat->copy_gates_from(oc_U, ids_U_target);

        // // --- Transfer the gates to the main circuit ---
        // auto all_qubits = YMATH::get_range(0, nq_);
        // copy_gates_from(
        //     oc_repeat,
        //     all_qubits,
        //     YSB(nullptr), 
        //     cs_unit, cs_zero,
        //     flag_inv        
        // ); 

        for(int i_mult = 0; i_mult < N_mult; i_mult++)
            copy_gates_from(
                oc_U_in,
                ids_U_target,
                YSB(nullptr), 
                cs_unit, cs_zero,
                flag_inv        
            ); 
    }
    timer.StopPrint();
    return get_the_circuit();
}





YQCP QCircuit::phase_estimation(
    YCVI ta, 
    const std::shared_ptr<const QCircuit>& A, 
    const std::shared_ptr<const QCircuit>& INIT,
    YCVI ty, 
    YCVI cs_unit, YCVI cs_zero,  
    YCB flag_inv,
    YCB flag_box
){
    string pe_name_tex = "PE";
    auto ny      = ty.size();
    auto nq_A    = A->get_n_qubits();
    if(ta.size() != nq_A)
    {
        string err_line;
        err_line  = "--- Error: setting the structure of the PE gate ---\n";
        err_line += "The subcircuit " + A->get_name() + " has " + to_string(nq_A) + " qubits, while" + 
            " one has indicated " + to_string(ta.size()) + " qubits to connect to.";
        throw err_line;
    }
    if(ta.size() != INIT->get_n_qubits())
    {
        string err_line;
        err_line  = "--- Error: setting the structure of the PE gate ---\n";
        err_line += "The subcircuit " + INIT->get_name() + " has " + to_string(INIT->get_n_qubits()) + " qubits, while" + 
            " one has indicated " + to_string(ta.size()) + " qubits to connect to.";
        throw err_line;
    }

    // --- create sequences from several A operators ---
    list<shared_ptr<QCircuit>> sequs_A;
    for(int iy = 0; iy < ny; iy++)
    {
        uint32_t N_rot = 1 << iy;
        auto c_temp = make_shared<QCircuit>(
            "A2^{"s+to_string(iy)+"}"s, env_, path_to_output_, nq_A
        );
        auto qa = c_temp->add_register("a", nq_A);
        for(int i_rot = 0; i_rot < N_rot; i_rot++)
            c_temp->copy_gates_from(A, qa);
        sequs_A.push_back(c_temp);
    }

    // --- create an envelop circuit for the phase estimation operator ---
    int count_c = -1;
    auto oc_pe = make_shared<QCircuit>(
        "PE", env_, path_to_output_, nq_A + ny
    );
    auto qy = oc_pe->add_register("y", ny);
    auto qa = oc_pe->add_register("a", nq_A);
    oc_pe->copy_gates_from(
        INIT, qa, 
        YMBo("INIT", qa)
        // YSB(nullptr)
    );
    oc_pe->h(qy);
    for(auto& one_sequ: sequs_A)
    {
        count_c++;
        auto ids_t_sequ = YVIv(qa);
        ids_t_sequ.push_back(qy[count_c]);
        oc_pe->copy_gates_from(
            one_sequ, 
            ids_t_sequ,  
            YMBo(one_sequ->get_name(), qa, YVIv{}, YVIv{}, one_sequ->get_name()),
            // YSB(nullptr), // without embedding into a box;
            {qy[count_c]}, YVIv {},
            false
        );
    }
    oc_pe->quantum_fourier(qy, YVIv{}, YVIv{}, true, true);

    // --- invert the circuit if necessary ---
    if(flag_inv)
    {
        oc_pe->h_adjoint();
        pe_name_tex += "^\\dagger";
    }

    // --- copy the PE env. circuit to the current circuit ---
    auto circ_qubits = YVIv({ta});
    circ_qubits.insert(circ_qubits.end(), ty.begin(), ty.end());

    auto box = YSB(nullptr);
    if(flag_box)
        box = YMBo("PE", circ_qubits, YVIv{}, YVIv{}, pe_name_tex);
    copy_gates_from(
        oc_pe,
        circ_qubits,
        box, 
        cs_unit, cs_zero,
        false      
    );
    return get_the_circuit();
}


qreal QCircuit::get_value_from_word(YCS word)
{
    if(word.find("<") == string::npos)
    {
        istringstream sstr(word);
        qreal res_value;
        if(!(sstr >> res_value))
            throw "Wrong format"s;
        return res_value;
    }

    unsigned first = word.find("<");
    unsigned last = word.find(">");
    string const_name = word.substr(first+1,last-first-1);

    if(constants_.find(const_name) == constants_.end())
        throw "The constant with the name "s + const_name + " is not found."s;
    return constants_[const_name];
}


void  QCircuit::qsvt_read_parameters(YCS gate_name, QSVT_pars& data)
{
    string filename = path_to_output_ + "/" + gate_name + FORMAT_QSP;
    ifstream ff_qsvt(filename);
    if(!ff_qsvt.is_open()) throw "Error: there is no file "s + filename;

    string line, key_name;
    while (getline(ff_qsvt, line))
    {
        key_name = "";
        line = YMIX::remove_comment(line);
        if(line.find_first_not_of(' ') == string::npos)
            continue;

        istringstream iss(line);
        iss >> key_name;

        if(YMIX::compare_strings(key_name, "filename_angles"))
        {
            iss >> data.filename_angles;
            continue;
        }
        if(YMIX::compare_strings(key_name, "n_repeat"))
        {
            iss >> data.n_repeat;
            continue;
        }
    }
    ff_qsvt.close();

    if(data.filename_angles.empty())
    {
        throw "Name of the file with QSVT angles is not indicated in "s + filename;
        exit(-1);
    }
        
    // --- read the .hdf5 file with angles ---
    string temp = data.filename_angles;
    if(!YMIX::compare_strings(
        temp.substr(
            temp.size()-string(FORMAT_HDF5).size(),string(FORMAT_HDF5).size()
        ), FORMAT_HDF5)
    ) temp += FORMAT_HDF5;
    data.filename_angles = temp;

    YMIX::print_log("\nRead angles from the file: "s + data.filename_angles);

    YMIX::H5File ff;
    ff.set_name(path_to_output_ + "/" + data.filename_angles);
    ff.open_r();

    // cout << endl;
    // if(ff.is_exist("basic")) cout << "Group basic exists." << endl;
    //     else cout << "Group basic does not exist." << endl;
    // if(not ff.is_exist("angles")) cout << "Group angles does not exist." << endl;
    //     else cout << "Group angles exists." << endl;

    // --- Old format of the .hdf5 file with angles ---
    if(ff.is_exist("angles"))
    {
        ff.read_scalar(data.type, "polynomial_type", "basic");
        ff.read_scalar(data.rescaling_factor, "rescaling_factor", "basic");
        ff.read_scalar(data.eps_qsvt, "eps", "basic");
        ff.read_scalar(data.f_par,    "par", "basic");
        if(YMIX::compare_strings(data.type, std::vector<std::string> {"matrix-inversion", "xgaussian"}))
        {
            ff.read_vector(data.angles_phis_odd, "odd", "angles");
            data.parity = 1;
        }
        else if(YMIX::compare_strings(data.type, "gaussian-arcsin"))
        {
            ff.read_vector(data.angles_phis_even, "even", "angles");
            data.parity = 0;
        }
        else if(YMIX::compare_strings(data.type, "QSVT-ham"))
        {
            ff.read_vector(data.angles_phis_odd,   "odd", "angles");
            ff.read_vector(data.angles_phis_even, "even", "angles");
            data.parity = -1;
        }
        else if(YMIX::compare_strings(data.type, "QSP-ham"))
        {
            // write all angles (odd and even) into the variable data.angles_phis_odd:
            ff.read_vector(data.angles_phis_arbitrary, "QSP-ham", "angles"); 
            data.parity = -1;
        }
        else
        {
            throw string("QSVT polynomial type " + data.type + " is not recognized.");
        }
    }
    // --- NEW format of the .hdf5 file with angles ---
    else
    {
        ff.read_scalar(data.type,   "function-type",      "basic");
        ff.read_scalar(data.f_par,  "function-parameter", "basic");
        ff.read_scalar(data.parity, "function-parity",    "basic");
        ff.read_scalar(data.rescaling_factor, "factor-norm", "basic");
        ff.read_scalar(data.eps_qsvt,         "abs-error",   "basic");
        if(data.parity == 0)
        {
            ff.read_vector(data.angles_phis_even, "phis", "results");
        }
        if(data.parity == 1)
        {
            ff.read_vector(data.angles_phis_odd, "phis", "results");
        }
    }
    ff.close();

    // --- print resulting parameters ---
    stringstream istr;
    istr << "--- QSVT gate with the name: " << gate_name << " ---\n";
    istr << "   QSVT type:  "            << data.type  << ";\n";
    istr << "   rescaling factor:  "     << data.rescaling_factor << ";\n";
    istr << "   QSVT error: "            << data.eps_qsvt << ";\n";
    istr << "   Polynomial parity: "     << data.parity << ";\n";
    istr << "   Polynomial parameter: "  << data.f_par << ";\n";
    if(data.parity == 0)
    {
        istr << "   number of angles: " << data.angles_phis_even.size() << ";\n";
    }
    if(data.parity == 1)
    {
        istr << "   number of angles: " << data.angles_phis_odd.size() << ";\n";
    }
    if(YMIX::compare_strings(data.type, "QSP-ham"))
    {
        istr << "   number of angles: " << data.angles_phis_arbitrary.size() << ";\n";
        istr << "   single time interval: " << data.f_par << ";\n";
        istr << "   number of time intervals: " << data.n_repeat << ";\n";
    }
    YMIX::print_log(istr.str());
}


YQCP QCircuit::qsvt_def_parity(
    YCVQ phis_in,
    YCI a_qsvt,
    YCVI qs_be_in, 
    const std::shared_ptr<const QCircuit> BE,
    YCVI cs_unit, YCVI cs_zero, 
    YCB flag_inv,
    YCB flag_box
){
    YMIX::YTimer timer;
    auto phis = YVQv(phis_in);
    auto N_angles = phis.size();
    auto n_be     = BE->get_n_qubits();
    auto n_be_anc = BE->get_na();
    auto cs_total_zero = YVIv(cs_zero);

    string qsvt_name_tex = "QSVT";
    string be_box_name = "BE";
    string be_box_name_tex;
    string be_box_cc_name_tex;

    // N_angles = 4; /// for testing;

    // --- separate BE ancillae and input qubits ---
    auto qs_be = YVIv(qs_be_in);
    vector<int> be_input(qs_be.begin(),                   qs_be.begin() + n_be - n_be_anc);
    vector<int>   be_anc(qs_be.begin() + n_be - n_be_anc, qs_be.end()                    );

    // --- pre-initialize the BE oracle ---
    // of the same size as the whole current circuit:
    auto oc_be = make_shared<QCircuit>(be_box_name, env_, path_to_output_, nq_);
    oc_be->add_register("r", nq_);
    oc_be->copy_gates_from(BE, qs_be, YSB(nullptr), cs_unit, cs_zero, flag_inv);

    // --- create the adjoint block-encoding oracle ---
    auto oc_be_inv = make_shared<QCircuit>(oc_be);
    oc_be_inv->h_adjoint();

    // --- QSVT circuit ---
    be_box_name_tex    = be_box_name;
    be_box_cc_name_tex = be_box_name + "^\\dagger"s;
    if(flag_inv)
    {
        reverse(phis.begin(), phis.end());
        be_box_name_tex    = be_box_cc_name_tex;
        be_box_cc_name_tex = be_box_name;
    }

    // to control on zero ancilla (and other zero-control nodes):
    cs_total_zero.insert(cs_total_zero.end(), be_anc.begin(), be_anc.end());
    sort(cs_total_zero.begin(), cs_total_zero.end());

    timer.StartPrint("Creating the QSVT circuit... ");
    h(a_qsvt, cs_unit, cs_zero);

    // form the first controlled projector (here, X-gates are necessary):
    x(a_qsvt, cs_unit, cs_total_zero);
    rz(a_qsvt, 2*phis[0], cs_unit, cs_zero, flag_inv);
    x(a_qsvt, cs_unit, cs_total_zero);

    // other controlled projectors:
    for(uint32_t count_angle = 1; count_angle < N_angles; ++count_angle)
    {
        // oracle:
        insert_gates_from(
            oc_be.get(), 
            // make_shared<Box__>(be_box_name, qs_be, cs_unit, cs_zero, be_box_name_tex)
            YSB(nullptr) // without embedding into a box
        );

        if((N_angles % 2) == 0)
        {
            if(!flag_inv)
            {
                if(count_angle == (N_angles - 1)) 
                    z(a_qsvt, cs_unit, cs_zero);
            }
            else
            {
                if(count_angle == 1) 
                    z(a_qsvt, cs_unit, cs_zero);
            }
        }

        // projector:
        x(a_qsvt, cs_unit, cs_total_zero);
        rz(a_qsvt, 2*phis[count_angle], cs_unit, cs_zero, flag_inv);
        x(a_qsvt, cs_unit, cs_total_zero);
        
        count_angle += 1;
        if(count_angle < N_angles)
        {
            // inverse oracle:
            insert_gates_from(
                oc_be_inv.get(), 
                // make_shared<Box__>(be_box_name, qs_be, cs_unit, cs_zero, be_box_cc_name_tex)
                YSB(nullptr) // without embedding into a box
            );

            // projector:
            x(a_qsvt, cs_unit, cs_total_zero);
            rz(a_qsvt, 2*phis[count_angle], cs_unit, cs_zero, flag_inv);
            x(a_qsvt, cs_unit, cs_total_zero);
        }
    }
    h(a_qsvt, cs_unit, cs_zero);
    timer.StopPrint();
    
    return get_the_circuit();
}


YQCP QCircuit::qetu_def_parity(
    YCVQ phis_in,
    YCI a_qsvt,
    YCVI qs_be_in, 
    const std::shared_ptr<const QCircuit> BE,
    YCVI cs_unit, YCVI cs_zero, 
    YCB flag_inv,
    YCB flag_box
){
    YMIX::YTimer timer;
    auto phis = YVQv(phis_in);
    auto N_angles = phis.size();
    auto n_be     = BE->get_n_qubits();
    auto n_be_anc = BE->get_na();

    string qsvt_name_tex = "QETU";
    string be_box_name = "BE";
    string be_box_name_tex;
    string be_box_cc_name_tex;

    // N_angles = 4; /// for testing;

    // --- separate BE ancillae and input qubits ---
    auto qs_be = YVIv(qs_be_in);
    vector<int> be_input(qs_be.begin(),                   qs_be.begin() + n_be - n_be_anc);
    vector<int>   be_anc(qs_be.begin() + n_be - n_be_anc, qs_be.end()                    );

    // 1-control qubits of the BE oracle:
    auto cs_unit_qsvt = YVIv(cs_unit);
    cs_unit_qsvt.push_back(a_qsvt);

    // --- pre-initialize the BE oracle ---
    // of the same size as the whole current circuit:
    auto oc_be = make_shared<QCircuit>(be_box_name, env_, path_to_output_, nq_);
    oc_be->add_register("r", nq_);
    oc_be->copy_gates_from(BE, qs_be, YSB(nullptr), cs_unit_qsvt, cs_zero, flag_inv);

    // --- create the adjoint block-encoding oracle ---
    auto oc_be_inv = make_shared<QCircuit>(oc_be);
    oc_be_inv->h_adjoint();

    // --- QSVT circuit ---
    be_box_name_tex    = be_box_name;
    be_box_cc_name_tex = be_box_name + "^\\dagger"s;
    if(flag_inv)
    {
        reverse(phis.begin(), phis.end());
        be_box_name_tex    = be_box_cc_name_tex;
        be_box_cc_name_tex = be_box_name;
    }

    timer.StartPrint("Creating the QETU circuit... ");
    // h(a_qsvt, cs_unit, cs_zero);  

    // form the first controlled projector (here, X-gates are necessary):
    rx(a_qsvt, -2*phis[0], cs_unit, cs_zero, flag_inv);

    // other controlled projectors:
    for(uint32_t count_angle = 1; count_angle < N_angles; ++count_angle)
    {
        // oracle:
        insert_gates_from(
            oc_be.get(), 
            // make_shared<Box__>(be_box_name, qs_be, cs_unit_qsvt, cs_zero, be_box_name_tex)
            YSB(nullptr) // without embedding into a box
        );
        rx(a_qsvt, -2*phis[count_angle], cs_unit, cs_zero, flag_inv);

        count_angle += 1;
        if(count_angle < N_angles)
        {
            // inverse oracle:
            insert_gates_from(
                oc_be_inv.get(), 
                // make_shared<Box__>(be_box_name, qs_be, cs_unit_qsvt, cs_zero, be_box_cc_name_tex)
                YSB(nullptr) // without embedding into a box
            );
            rx(a_qsvt, -2*phis[count_angle], cs_unit, cs_zero, flag_inv);
        }
    }
    // h(a_qsvt, cs_unit, cs_zero);
    timer.StopPrint();
    
    return get_the_circuit();
}



YQCP QCircuit::qsp_ham(
        YCS name_qsp_circuit,
        YCVQ angles_qsp,
        YCI nt,
        YCI a_qsp,
        YCI a_qu, 
        YCVI qs_be_in, 
        const std::shared_ptr<const QCircuit> BE,
        YCVI cs_unit, 
        YCVI cs_zero, 
        YCB flag_inv
){
    YMIX::YTimer timer;
    auto phis = YVQv(angles_qsp);
    auto N_angles = phis.size();
    auto n_be     = BE->get_n_qubits();
    auto n_be_anc = BE->get_na();
    auto all_qubits = YMATH::get_range(0, nq_);

    // N_angles = 3; // for testing;

    timer.StartPrint("Creating the QSP circuit (OPT1)... ");

    // --- separate BE ancillae and input qubits ---
    auto qs_be = YVIv(qs_be_in);
    YVIv be_input(qs_be.begin(),                 qs_be.begin() + n_be - n_be_anc);
    YVIv be_anc(qs_be.begin() + n_be - n_be_anc, qs_be.end()                    );

    // --- qubitization oracle W ---
    auto oW = make_shared<QCircuit>("W", env_, path_to_output_, nq_);
    //     std::map<std::string, qreal>(), false, true 
    // );
    oW->add_register("r", nq_);

    // block-encoding oracles:
    oW->copy_gates_from(BE, qs_be, YSB(nullptr), YVIv{},     YVIv{a_qu});
    oW->copy_gates_from(BE, qs_be, YSB(nullptr), YVIv{a_qu}, YVIv{}, true);

    // reflector:
    oW->h(a_qu);
    oW->phase_zero(a_qu, M_PI, YVIv{}, be_anc);
    oW->phase_zero(a_qu, M_PI);
    oW->h(a_qu);

    // oW->save_regs();
    // oW->print_gates();

    // --- QSP circuit for one interval ---
    qreal aa; 
    auto qubits_W = YVIv(qs_be);
    qubits_W.push_back(a_qu);

    auto oc_qsp = make_shared<QCircuit>("QSP", env_, path_to_output_, nq_);
    //     std::map<std::string, qreal>(), false, true 
    // );
    oc_qsp->add_register("r", nq_);

    oc_qsp->h(a_qsp)->h(a_qu);
    for(unsigned count_angle = 0; count_angle < int((N_angles-1)/2); ++count_angle)
    {
        aa = phis[2*count_angle];
        oc_qsp->rz(a_qsp, -aa)->h(a_qsp);
        oc_qsp->copy_gates_from(oW, all_qubits, YSB(nullptr), YVIv{}, YVIv{a_qsp}); 
        oc_qsp->phase_zero(a_qsp, -M_PI_2);
        oc_qsp->h(a_qsp)->rz(a_qsp, aa);

        aa = phis[2*count_angle+1];
        oc_qsp->rz(a_qsp, -aa)->h(a_qsp);
        oc_qsp->copy_gates_from(oW, all_qubits, YSB(nullptr), YVIv{a_qsp}, YVIv{}, true); 
        oc_qsp->phase(a_qsp, M_PI_2);
        oc_qsp->h(a_qsp)->rz(a_qsp, aa);
    }
    aa = phis[N_angles-1];
    oc_qsp->rz(a_qsp, aa);
    oc_qsp->h(a_qsp)->h(a_qu);


    // oc_qsp->save_regs();
    // oc_qsp->print_gates();

    // --- Several time intervals ---
    auto qsp_qubits = YVIv(qubits_W);
    qsp_qubits.push_back(a_qsp);

    // to save the state just before the QSP simulation:
    string name_stop;
    if(flag_stop_gates_)
    {
        name_stop = "QSP-H <" + name_qsp_circuit + ">: t = " + to_string(0);
        add_stop_gate(name_stop);
        // cout << "QSP: add stop gate" << endl;
    }
        
    // QSP simulation with several intervals
    for(int it = 0; it < nt; it++)
    { 
        copy_gates_from(
            oc_qsp,
            all_qubits,
            YSB(nullptr), 
            cs_unit, cs_zero,
            flag_inv      
        ); 

        // to save a state after the simulation of a time interval:
        if(flag_stop_gates_)
        {
            name_stop = "QSP-H <" + name_qsp_circuit + ">: t = " + to_string(it+1);
            add_stop_gate(name_stop);
            // cout << "QSP: add stop gate" << endl;
        }
    }
    timer.StopPrint();
    return get_the_circuit();
}



YQCP QCircuit::qsp_ham_opt2(
        YCS name_qsp_circuit,
        YCVQ angles_qsp,
        YCI nt,
        YCI a_qsp,
        YCI a_qu, 
        YCVI qs_be_in, 
        const std::shared_ptr<const QCircuit> BE,
        YCVI cs_unit, 
        YCVI cs_zero, 
        YCB flag_inv
){
    YMIX::YTimer timer;
    auto phis = YVQv(angles_qsp);
    auto N_angles = phis.size();
    auto n_be     = BE->get_n_qubits();
    auto n_be_anc = BE->get_na();
    auto all_qubits = YMATH::get_range(0, nq_);

    timer.StartPrint("Creating the QSP circuit (OPT2)... ");

    // --- separate BE ancillae and input qubits ---
    auto qs_be = YVIv(qs_be_in);
    YVIv be_input(qs_be.begin(),                 qs_be.begin() + n_be - n_be_anc);
    YVIv be_anc(qs_be.begin() + n_be - n_be_anc, qs_be.end()                    );

    // --- qubitization oracle W ---
    auto oW = make_shared<QCircuit>("W", env_, path_to_output_, nq_);

    // block-encoding oracles:
    oW->copy_gates_from(BE, qs_be, YSB(nullptr), YVIv{},     YVIv{a_qu});
    oW->copy_gates_from(BE, qs_be, YSB(nullptr), YVIv{a_qu}, YVIv{}, true);

    // reflector:
    oW->h(a_qu);
    oW->phase_zero(a_qu, M_PI, YVIv{}, be_anc);
    oW->phase_zero(a_qu, M_PI);
    oW->h(a_qu);

    // inverse oW
    auto oWi = make_shared<QCircuit>(oW);
    oWi->h_adjoint();

    // add control qubits
    auto ctrl_0_ow = YVIv(cs_zero);
    ctrl_0_ow.push_back(a_qsp);
    oW->add_control_qubits(cs_unit, ctrl_0_ow);

    auto ctrl_1_owi = YVIv(cs_unit);
    ctrl_1_owi.push_back(a_qsp);
    oWi->add_control_qubits(ctrl_1_owi, cs_zero);

    // to save the state just before the QSP simulation:
    string name_stop;
    if(flag_stop_gates_)
    {
        name_stop = "QSP-H <" + name_qsp_circuit + ">: t = " + to_string(0);
        add_stop_gate(name_stop);
    }
        
    // QSP simulation with several intervals
    qreal aa; 
    for(int it = 0; it < nt; it++)
    { 
        if(flag_inv)
        {
            h(a_qsp, cs_unit, cs_zero)->h(a_qu, cs_unit, cs_zero);
            aa = phis[N_angles-1];
            rz(a_qsp, -aa, cs_unit, cs_zero);
            for(int count_angle = int((N_angles-1)/2)-1; count_angle >= 0; --count_angle)
            {
                aa = phis[2*count_angle+1];
                rz(a_qsp, -aa, cs_unit, cs_zero);
                h(a_qsp, cs_unit, cs_zero);
                phase(a_qsp, -M_PI_2, cs_unit, cs_zero);
                insert_gates_from(oW.get());
                h(a_qsp, cs_unit, cs_zero);
                rz(a_qsp, aa, cs_unit, cs_zero);

                aa = phis[2*count_angle];
                rz(a_qsp, -aa, cs_unit, cs_zero);
                h(a_qsp, cs_unit, cs_zero);
                phase_zero(a_qsp, M_PI_2, cs_unit, cs_zero);
                insert_gates_from(oWi.get());
                h(a_qsp, cs_unit, cs_zero);
                rz(a_qsp, aa, cs_unit, cs_zero);
            }
            h(a_qsp, cs_unit, cs_zero)->h(a_qu, cs_unit, cs_zero);
        }
        else
        {
            h(a_qsp, cs_unit, cs_zero)->h(a_qu, cs_unit, cs_zero);
            for(unsigned count_angle = 0; count_angle < int((N_angles-1)/2); ++count_angle)
            {
                aa = phis[2*count_angle];
                rz(a_qsp, -aa, cs_unit, cs_zero);
                h(a_qsp, cs_unit, cs_zero);
                insert_gates_from(oW.get());
                phase_zero(a_qsp, -M_PI_2, cs_unit, cs_zero);
                h(a_qsp, cs_unit, cs_zero);
                rz(a_qsp, aa, cs_unit, cs_zero);

                aa = phis[2*count_angle+1];
                rz(a_qsp, -aa, cs_unit, cs_zero);
                h(a_qsp, cs_unit, cs_zero);
                insert_gates_from(oWi.get());
                phase(a_qsp, M_PI_2, cs_unit, cs_zero);
                h(a_qsp, cs_unit, cs_zero);
                rz(a_qsp, aa, cs_unit, cs_zero);
            }
            aa = phis[N_angles-1];
            rz(a_qsp, aa, cs_unit, cs_zero);
            h(a_qsp, cs_unit, cs_zero)->h(a_qu, cs_unit, cs_zero);
        }

        // to save a state after the simulation of a time interval:
        if(flag_stop_gates_)
        {
            name_stop = "QSP-H <" + name_qsp_circuit + ">: t = " + to_string(it+1);
            add_stop_gate(name_stop);
        }
    }
    timer.StopPrint();
    return get_the_circuit();
}



YQCP QCircuit::selector_power(
    YCVI rs, 
    const std::shared_ptr<const QCircuit> oc_U, YCVI ids_U_target,
    YCVI cs_unit, YCVI cs_zero,
    YCB flag_inv
){
    YMIX::YTimer timer;
    timer.StartPrint("Creating the SelectorPower gate... ");

    uint32_t ns, is_res;
    ns = rs.size();
    for(uint32_t is = 0; is < ns; is++)
    {
        if(flag_inv)
            is_res = ns - is - 1;
        else
            is_res = is;
        auto loc_cs_unit = YVIv(cs_unit);
        loc_cs_unit.push_back(rs[is_res]);
        for(uint32_t i_repeat = 0; i_repeat < (1<<is_res); i_repeat++)
            copy_gates_from(oc_U, ids_U_target, YSB(nullptr), loc_cs_unit, cs_zero, flag_inv);
    }
    timer.StopPrint();
    return get_the_circuit();
}



// YQCP QCircuit::selector_power(
//     YCVI rs, 
//     const std::shared_ptr<const QCircuit> oc_U, YCVI ids_U_target,
//     YCVI cs_unit, YCVI cs_zero,
//     YCB flag_inv
// ){
//     YMIX::YTimer timer;
//     timer.StartPrint("Creating the SelectorPower gate... ");
//     std::shared_ptr<QCircuit> oc_temp;
//     uint32_t ns, Ns;
//     auto all_qubits = YMATH::get_range(0, nq_);

//     ns = rs.size();
//     Ns = 1 << ns;
//     auto oc_selector = make_shared<QCircuit>("S", env_, path_to_output_, nq_);

//     for(uint32_t is = 0; is < ns; is++)
//         for(uint32_t i_repeat = 0; i_repeat < (1<<is); i_repeat++)
//             oc_selector->copy_gates_from(oc_U, ids_U_target, YSB(nullptr), YVIv{rs[is]});

//     // --- Transfer the gates to the main circuit ---
//     copy_gates_from(
//         oc_selector,
//         all_qubits,
//         YSB(nullptr), 
//         cs_unit, cs_zero,
//         flag_inv        
//     ); 
//     timer.StopPrint();
//     return get_the_circuit();
// }



YQCP QCircuit::LCHS_QSP(
    const std::shared_ptr<const QCircuit> Ut, YCVI ids_Ut, YCI Nt,
    const std::shared_ptr<const QCircuit> Ow, YCVI ids_Ow,
    YCVI cs_unit, YCVI cs_zero,
    YCB flag_inv
){
    auto all_qubits = YMATH::get_range(0, nq_);
    
    // --- get ancilla qubits of the oracles Ut and Ow ---
    auto qs_Ut = YVIv(ids_Ut);
    int n_anc_Ut = Ut->get_na();
    int n_inp_Ut = Ut->get_n_qubits() - n_anc_Ut;

    auto qs_Ow = YVIv(ids_Ow);
    int n_anc_Ow = Ow->get_na();
    int n_inp_Ow = Ow->get_n_qubits() - n_anc_Ow;

    YVIv ids_anc_Ow(qs_Ow.begin()+n_inp_Ow, qs_Ow.end());
    YVIv ids_anc_Ut(qs_Ut.begin()+n_inp_Ut, qs_Ut.end());

    // - ancilla qubits of Ut that do not coincide with the ancillae of Ow -
    YVIv ids_anc_Ut_unique(ids_anc_Ut.size());
    int temp;
    int counter = -1;
    for(int ii = 0; ii < ids_anc_Ut.size(); ii++)
    {
        temp = ids_anc_Ut[ii];
        if(
            find(ids_Ow.begin(), ids_Ow.end(), temp) == ids_Ow.end()
        ){
            counter++;
            ids_anc_Ut_unique[counter] = temp;
        }
    }
    ids_anc_Ut_unique.resize(counter+1);

    // --- Create the LCU selector (without using a compression gadget) ---
    auto oc_selector = make_shared<QCircuit>("selector", env_, path_to_output_, nq_);
    auto oc_temp_Ut = make_shared<QCircuit>("temp", env_, path_to_output_, nq_);
    oc_temp_Ut->copy_gates_from(Ut, ids_Ut);
    for(int it = 0; it < Nt; it++)
        oc_selector->insert_gates_from(oc_temp_Ut.get());

    // --- Data for the compression gadget ---
    GADGET_pars data_cg;

    // --- Create the LCHS (LCU) circuit ---
    auto oc_LCHS = make_shared<QCircuit>("LCHS-QSP", env_, path_to_output_, nq_);

    // - left Ow -
    oc_LCHS->copy_gates_from(Ow,          ids_Ow,     YSB(nullptr), YVIv{}, ids_anc_Ut_unique, false);

    // - add the selector -
    oc_LCHS->copy_gates_from(oc_selector, all_qubits, YSB(nullptr), YVIv{}, ids_anc_Ow);

    // // - selector via the compression gadget -
    // compression_gadget(
    //     data_cg, ids_counter, oc_U, ids_U_target, N_mult, flag_step_output, 
    //     ids_unit, ids_zero, flag_inv
    // );

    // - right Ow-adjoint -
    oc_LCHS->copy_gates_from(Ow,          ids_Ow,     YSB(nullptr), YVIv{}, ids_anc_Ut_unique, true);

    // --- Transfer the gates to the main circuit ---
    copy_gates_from(
        oc_LCHS,
        all_qubits,
        YSB(nullptr), 
        cs_unit, cs_zero,
        flag_inv        
    ); 
    return get_the_circuit();
}



YQCP QCircuit::DirDec_Y(
        YCVQ phis_y_in,
        YCVI ids_ctrl_in, 
        YCI id_targ,
        YCVI cs_unit, YCVI cs_zero,
        YCB flag_inv
){
    YVIv ids_ctrl = vector<int>(ids_ctrl_in);
    YVQv phis_y   = vector<qreal>(phis_y_in);
    int nc = ids_ctrl.size();
    int Nc = 1 << nc;

    if(flag_inv)
        std::reverse(phis_y.begin(), phis_y.end());
    std::sort(ids_ctrl.begin(), ids_ctrl.end());
    for(int i_int = 0; i_int < Nc; i_int++)
    {
        YVshv bs(nc);

        int i_int_new = i_int;
        if(flag_inv)
            i_int_new = Nc - i_int - 1;
        YMATH::intToBinary(i_int_new, bs);
        // cout << i_int << endl;
        // YMIX::print(bs);

        YVIv ids_unit = vector<int>(cs_unit);
        YVIv ids_zero = vector<int>(cs_zero);
        for(int ib = 0; ib < nc; ib++)
        {
            if(bs[ib] == 0)
                ids_zero.push_back(ids_ctrl[nc - ib - 1]);
            else
                ids_unit.push_back(ids_ctrl[nc - ib - 1]);
        }
        // YMIX::print(ids_zero);
        // YMIX::print(ids_unit);
        // cout << endl;

        ry(id_targ, phis_y[i_int], ids_unit, ids_zero, flag_inv);
    }
    return get_the_circuit();
}


YQCP QCircuit::DirDec_C(
        YCVQ phis_y_in,
        YCVQ phis_z_in,
        YCVI ids_ctrl_in, 
        YCI id_targ,
        YCVI cs_unit, YCVI cs_zero,
        YCB flag_inv
){
    YVIv ids_ctrl = vector<int>(ids_ctrl_in);
    YVQv phis_y   = vector<qreal>(phis_y_in);
    YVQv phis_z   = vector<qreal>(phis_z_in);
    int nc = ids_ctrl.size();
    int Nc = 1 << nc;
    if(flag_inv)
    {
        std::reverse(phis_y.begin(), phis_y.end());
        std::reverse(phis_z.begin(), phis_z.end());
    }
        
    std::sort(ids_ctrl.begin(), ids_ctrl.end());
    for(int i_int = 0; i_int < Nc; i_int++)
    {
        YVshv bs(nc);
        
        int i_int_new = i_int;
        if(flag_inv)
            i_int_new = Nc - i_int - 1;
        YMATH::intToBinary(i_int_new, bs);

        YVIv ids_unit = vector<int>(cs_unit);
        YVIv ids_zero = vector<int>(cs_zero);
        for(int ib = 0; ib < nc; ib++)
        {
            if(bs[ib] == 0)
                ids_zero.push_back(ids_ctrl[nc - ib - 1]);
            else
                ids_unit.push_back(ids_ctrl[nc - ib - 1]);
        }
        rc(id_targ, phis_z[i_int], phis_y[i_int], ids_unit, ids_zero, flag_inv);
    }
    return get_the_circuit();
}





void QCircuit::get_nonancilla_regs(
    std::vector<std::string>& reg_names_nonanc,
    YVI n_qubits
){
    reg_names_nonanc.clear();
    n_qubits.clear();
    for(int ireg = 0; ireg < regnames_.size(); ireg++)
    {
        string reg_name = regnames_[ireg];
        if(!flags_anc_regs_[reg_name]) // if nonancilla, then...
        {
            reg_names_nonanc.push_back(reg_name);
            n_qubits.push_back(regs_[reg_name].size());
        }
    }
}

