#ifndef QGATES_H
#define QGATES_H

#include "QLib.h"
#include "QuEST_lib.h"


// -----------------------------------------------------------------
// --- Gates ---
// -----------------------------------------------------------------

class Gate__
{
    public:
        Gate__(YCS name) : name_(name)
        {
            id_layer_ = -1;
            flag_h_adjoint_ = false;
            tex_name_ = name_;
        }
        ~Gate__(){}

        /**
         * @brief Copy Construct.
         */
        Gate__(const Gate__& oo)
        {
            // copy to this object:

            name_ = oo.name_;
            type_ = oo.type_;
            tex_name_ = oo.tex_name_;
            
            ts_ = YVIv(oo.ts_);
            cs_unit_ = YVIv(oo.cs_unit_);
            cs_zero_ = YVIv(oo.cs_zero_);

            copy_matrix2(oo.u2_, this->u2_);
            if(oo.un_a_ != nullptr)
                un_a_ = std::make_shared<YMATH::YMatrix>(oo.un_a_);
            if(oo.un_b_ != nullptr)
                un_b_ = std::make_shared<YMATH::YMatrix>(oo.un_b_);

            pars_ = YVQv(oo.pars_);

            add_inf_ = std::string(oo.add_inf_);

            id_layer_ = oo.get_layer();

            flag_h_adjoint_ = oo.flag_h_adjoint_;
            flag_start_ = oo.flag_start_;
        }

        // copy this gate to a new gate:
        virtual YSG copy_gate() const { return std::make_shared<Gate__>(*this); }
        virtual void h_adjoint(){ flag_h_adjoint_ = !flag_h_adjoint_; };
        virtual void generate(Qureg& oc){};
        virtual void write_to_file(YMIX::File& cf){ write_to_file_base(cf); };

        virtual void activate_gadget(YCI id_in_counter, YCI N_mult){};


        inline void check_control_nodes(YCVI control_qubits)
        {
            for(auto const& new_c: control_qubits)
            {
                auto it_unit = find(cs_unit_.begin(), cs_unit_.end(), new_c);
                auto it_zero = find(cs_zero_.begin(), cs_zero_.end(), new_c);
                if (it_unit != cs_unit_.end())
                {
                    std::ostringstream sstr;
                    sstr << "The gate " << name_ << " is already unit-controlled by the qubit [" << new_c << "]";
                    throw sstr.str();
                }
                if (it_zero != cs_zero_.end())
                {
                    std::ostringstream sstr;
                    sstr << "The gate " << name_ << " is already zero-controlled by the qubit [" << new_c << "]";
                    throw sstr.str();
                }
            } 
        }


        /**
         * @brief Add control qubits to the gate.
         * @param[in] unit_control_qubits 1-control nodes to add.
         * @param[in] zero_control_qubits 0-control nodes to add.
         */
        inline void add_control_qubits(YCVI unit_control_qubits, YCVI zero_control_qubits = {})
        {
            // unit-control nodes:
            check_control_nodes(unit_control_qubits);
            copy(unit_control_qubits.begin(), unit_control_qubits.end(), back_inserter(cs_unit_));

            // zero-control nodes:
            check_control_nodes(zero_control_qubits);
            copy(zero_control_qubits.begin(), zero_control_qubits.end(), back_inserter(cs_zero_));
        }

        inline std::string get_name(){ return name_; }

        /**
         * @brief new_pos[id-qubit-position] = id-qubit-new-position
         * E.g., assume a circuit has 3 qubits with indices 0,1,2.
         * If @param new_pos = [4,8,6], then new qubit positions are:
         * 0-> 4, 1->8, 2->6 (only possible if the niew circuit has al least 8+1 qubits).
         * @param[in] new_pos new positions of qubits. 
         */
        void correct_qubits(YCVI new_pos);

        inline 
        std::string get_type(){ return type_; }

        void get_gubits_act_on(YVI ids_qubit_act_on)
        {
            ids_qubit_act_on = YVIv(ts_);
            if(!cs_unit_.empty())
                ids_qubit_act_on.insert(ids_qubit_act_on.end(), cs_unit_.begin(), cs_unit_.end());
            if(!cs_zero_.empty())
                ids_qubit_act_on.insert(ids_qubit_act_on.end(), cs_zero_.begin(), cs_zero_.end());
        }

        void set_layer(const int64_t& id_layer){ id_layer_ = id_layer; }

        int64_t get_layer() const { return id_layer_; }

        inline void set_flag_start(YCB flag_start){ flag_start_ = flag_start; }
        inline bool get_flag_start(){ return flag_start_; }
        inline bool get_flag_h_adjoint(){ return flag_h_adjoint_; }

        void get_target_qubits(YVI ids_t){ ids_t = YVIv(ts_); }
        void get_unit_control_qubits(YVI ids_c){ ids_c = YVIv(cs_unit_); }
        void get_zero_control_qubits(YVI ids_c){ ids_c = YVIv(cs_zero_); }

        virtual void write_tex(
            std::vector<std::vector<std::string>>& tex_lines, 
            const int64_t& id_layer,
            YCU nq
        );

    protected:
        inline
        void change_element(YVI v_new, YCVI v, YCI old_x, YCI new_x)
        {
            int index;
            auto it = find(v.begin(), v.end(), old_x);
            if (it != v.end())
            {
                index = it - v.begin();
                v_new[index] = new_x;        
            }    
        }

        /**
         * Set a multi-controlled unitary in the QuEST circuit.
        */
        inline
        void mc_st_u(Qureg& oc, YCI t, YVI cs_unit, YVI cs_zero, const ComplexMatrix2& u)
        {
            multiMixControlledUnitary(oc, &cs_unit[0], cs_unit.size(), t, u, &cs_zero[0], cs_zero.size());
        }

        inline
        void copy_matrix2(const ComplexMatrix2& copy_from, ComplexMatrix2& copy_to)
        {
            copy_to.real[0][0] = copy_from.real[0][0]; 
            copy_to.real[0][1] = copy_from.real[0][1];
            copy_to.real[1][0] = copy_from.real[1][0];
            copy_to.real[1][1] = copy_from.real[1][1];

            copy_to.imag[0][0] = copy_from.imag[0][0]; 
            copy_to.imag[0][1] = copy_from.imag[0][1];
            copy_to.imag[1][0] = copy_from.imag[1][0];
            copy_to.imag[1][1] = copy_from.imag[1][1];
        }

        inline
        void write_to_file_base(YMIX::File& cf, const bool& flag_new_line=true)
        {
            cf << name_ << " ";
            
            if(flag_h_adjoint_)
                cf << "conj ";
            else
                cf << "orig "; 

            cf << "layer " << id_layer_ << " ";
            
            cf << "targets " << ts_.size() << " ";
            for(auto& x: ts_)
                cf << x << " ";

            cf << "controls " << cs_unit_.size() << " ";
            for(auto& x: cs_unit_)
                cf << x << " ";

            cf << "ocontrols " << cs_zero_.size() << " ";
            for(auto& x: cs_zero_)
                cf << x << " ";

            cf << "pars " << pars_.size() << " ";
            for(auto& x: pars_)
                cf << std::scientific << std::setprecision(PARAMETER_ACCURACY) << x << " ";
            
            if(flag_new_line)
                cf << "\n";
        }

        inline
        std::string tex_gate_width(
            std::vector<std::vector<std::string>>& tex_lines, 
            const int64_t& id_layer,
            YCU nq,
            YCI id_top_q
        ){
            std::string l_nq_gate;
            if(ts_.size() == 1)
                l_nq_gate = "";
            else
                l_nq_gate = "[" + std::to_string(abs(ts_.back() - ts_[0]+1)) + "]";
            return l_nq_gate;
        }

        inline
        std::string tex_get_gate_name(
            std::vector<std::vector<std::string>>& tex_lines, 
            const int64_t& id_layer,
            YCU nq,
            YCS l_brackets = "()",
            YCB flag_inv_par = false
        ){
            std::string l_name, l_par;
            l_name = tex_name_;
            if(!flag_inv_par && flag_h_adjoint_)
                l_name += std::string("^\\dagger");

            if(pars_.size())
            {
                l_name.push_back(l_brackets[0]);
                for(auto id_p = 0; id_p < pars_.size(); id_p++)
                {
                    std::stringstream sstr;
                    if(flag_inv_par && flag_h_adjoint_)
                        sstr << -pars_[id_p];
                    else
                        sstr << pars_[id_p];
                    sstr >> l_par;
                    if(id_p > 0)
                        l_name += std::string(", ");
                    l_name += l_par;
                }
                l_name.push_back(l_brackets[1]);
            }
            return l_name;
        }

        inline
        void tex_add_control(
            std::vector<std::vector<std::string>>& tex_lines, 
            const int64_t& id_layer,
            YCU nq,
            YCI id_top_q
        ){
            std::string l_c_dir;
            for(auto const& id_cq: cs_unit_)
            {
                l_c_dir = std::to_string(id_cq - id_top_q);
                tex_lines[nq - id_cq - 1][id_layer] = "&\\ctrl{" + l_c_dir + "}";
            }
            for(auto const& id_cq: cs_zero_)
            {
                l_c_dir = std::to_string(id_cq - id_top_q);
                tex_lines[nq - id_cq - 1][id_layer] = "&\\octrl{" + l_c_dir + "}";
            }
        }

        inline
        int get_most_signif_target_qubit(){ return *(std::max_element(ts_.begin(), ts_.end())); }

    public:
        const static std::string name_shared_;
        static const std::vector<std::string> avail_gate_names_;
        
    protected:
        std::string name_;// current name of the gate;
        std::string type_;// type of the gate;
        std::string tex_name_;// gate name as it is shown in the. tex file;

        YVIv ts_; // target qubits;
        YVIv cs_unit_; // control qubits;
        YVIv cs_zero_; // zero-control qubits;

        ComplexMatrix2 u2_; // matrix for a single-target gate;
        YVQv pars_; // parameters of the gate (e.g. angles);

        std::string add_inf_; // string line with additional information;

        YSM un_a_ = nullptr;
        YSM un_b_ = nullptr;

        int64_t id_layer_; // id of the circuit layer, where the gate sits on;

        bool flag_h_adjoint_; // whether the gate is Hermitiain adjoint or not;
        bool flag_start_ = true; // is it the left side of the box? 
};


class GStop__ : public Gate__
{
public:
    GStop__(YCS name) : Gate__(name){ type_ = "stop"; }
    YSG copy_gate() const { return std::make_shared<GStop__>(*this); }
    void write_to_file(YMIX::File& cf){}
    void write_tex(
            std::vector<std::vector<std::string>>& tex_lines, 
            const int64_t& id_layer,
            YCU nq
    ){}
};


class Box__ : public Gate__
{
    public:
        Box__(YCS name, YCVI ts, YCVI cs_unit = YVIv{}, YCVI cs_zero = YVIv{}, YCS tex_name="") : Gate__(name)
        {
            type_ = "box";
            ts_ = ts;
            cs_unit_ = cs_unit;
            cs_zero_ = cs_zero;

            if(!tex_name.empty())
                tex_name_ = tex_name;
        }

        YSG copy_gate() const {return std::make_shared<Box__>(*this);}
        YSB copy_box() const {return std::make_shared<Box__>(*this);}

        void write_to_file(YMIX::File& cf) override
        {
            cf << "Box ";
            if(flag_start_) 
                cf << "start ";
            else 
                cf << "end ";
            write_to_file_base(cf);
        }
};


class SQGate__ : public Gate__
{
    public:
        SQGate__(YCS name, YCI t) : Gate__(name)
        {   
            type_ = "q1";
            ts_ = YVIv(1);
            ts_[0] = t;
        }

        YSG copy_gate() const {return std::make_shared<SQGate__>(*this);};

        void h_adjoint(){
            flag_h_adjoint_ = !flag_h_adjoint_;
            u2_ = YMATH::inv_matrix2(u2_);
        }

        void generate(Qureg& oc){}
};


class X__ : public SQGate__
{
    public:
        X__(YCI t) : SQGate__(name_shared_, t){ u2_ = YGV::mX; }

        YSG copy_gate() const { return std::make_shared<X__>(*this); };

        void generate(Qureg& oc) 
        { 
            if(cs_unit_.empty() && cs_zero_.empty())
                pauliX(oc, ts_[0]);
            else
                mc_st_u(oc, ts_[0], cs_unit_, cs_zero_, u2_); 
        }

        void h_adjoint(){}

        void write_tex(
            std::vector<std::vector<std::string>>& tex_lines, 
            const int64_t& id_layer,
            YCU nq
        ){
            std::string l_nq_gate;

            // the most-signficant target qubit:
            auto id_top_q = get_most_signif_target_qubit();

            // gate keyword:
            if(cs_unit_.size() > 0 || cs_zero_.size() > 0)
                tex_lines[nq - id_top_q - 1][id_layer] = "&\\targ{}";
            else
                tex_lines[nq - id_top_q - 1][id_layer] = "&\\gate{X}";

            // add the control qubits 
            tex_add_control(tex_lines, id_layer, nq, id_top_q);
        }

    public:
        const static std::string name_shared_;
};


class Y__ : public SQGate__
{
    public:
        Y__(YCI t) : SQGate__(name_shared_, t){ u2_ = YGV::mY; }

        YSG copy_gate() const { return std::make_shared<Y__>(*this); };

        void generate(Qureg& oc) 
        { 
            if(cs_unit_.empty() && cs_zero_.empty())
                pauliY(oc, ts_[0]);
            else
                mc_st_u(oc, ts_[0], cs_unit_, cs_zero_, u2_); 
        }

        void h_adjoint(){}

    public:
        const static std::string name_shared_;
};


class Z__ : public SQGate__
{
    public:
        Z__(YCI t) : SQGate__(name_shared_, t){ u2_ = YGV::mZ; }

        YSG copy_gate() const { return std::make_shared<Z__>(*this); };

        void generate(Qureg& oc) 
        { 
            if(cs_unit_.empty() && cs_zero_.empty())
                pauliZ(oc, ts_[0]);
            else
                mc_st_u(oc, ts_[0], cs_unit_, cs_zero_, u2_); 
        }

        void h_adjoint(){}
    
    public:
        const static std::string name_shared_;
};


class H__ : public SQGate__
{
    public:
        H__(YCI t) : SQGate__(name_shared_, t){ u2_ = YGV::mH; }

        YSG copy_gate() const { return std::make_shared<H__>(*this); };

        void generate(Qureg& oc) 
        { 
            if(cs_unit_.empty() && cs_zero_.empty())
                hadamard(oc, ts_[0]); 
            else 
                mc_st_u(oc, ts_[0], cs_unit_, cs_zero_, u2_);
        }

        void h_adjoint(){}

    public:
        const static std::string name_shared_;
};


class sR__ : public SQGate__ // single-angle rotation
{
    public:
        sR__(YCS name, YCI t, YCQR a) : SQGate__(name, t){ pars_.push_back(a); }
        YSG copy_gate() const { return std::make_shared<sR__>(*this); };
        void generate(Qureg& oc) 
        { 
            if(cs_unit_.empty() && cs_zero_.empty())
                unitary(oc, ts_[0], u2_); // take a general unitary function since the u2_ can be inversed;
            else 
                mc_st_u(oc, ts_[0], cs_unit_, cs_zero_, u2_);
        }
};


class Rx__ : public sR__
{
    public:
        Rx__(YCI t, YCQR a) : sR__(name_shared_, t, a){ u2_ = YGV::mRx(a); tex_name_ = "R_x"; }
        YSG copy_gate() const { return std::make_shared<Rx__>(*this); };

    public:
        const static std::string name_shared_;
};


class Ry__ : public sR__
{
    public:
        Ry__(YCI t, YCQR a) : sR__(name_shared_, t, a){ u2_ = YGV::mRy(a); tex_name_ = "R_y"; }
        YSG copy_gate() const { return std::make_shared<Ry__>(*this); };

    public:
        const static std::string name_shared_;
};


class Rz__ : public sR__
{
    public:
        Rz__(YCI t, YCQR a) : sR__(name_shared_, t, a){ u2_ = YGV::mRz(a); tex_name_ = "R_z";}
        YSG copy_gate() const { return std::make_shared<Rz__>(*this); };

    public:
        const static std::string name_shared_;
};


class Phase__ : public sR__
{
    public:
        Phase__(YCI t, YCQR a) : sR__(name_shared_, t, a){ u2_ = YGV::mPhase(a); tex_name_ = "P"; }
        YSG copy_gate() const { return std::make_shared<Phase__>(*this); };

    public:
        const static std::string name_shared_;
};


class PhaseZero__ : public sR__
{
    public:
        PhaseZero__(YCI t, YCQR a) : sR__(name_shared_, t, a)
        { 
            u2_ = YGV::mPhaseZero(a); 
            tex_name_ = "P_0"; 
        }
        YSG copy_gate() const { return std::make_shared<PhaseZero__>(*this); };

    public:
        const static std::string name_shared_;
};


class sR2__ : public sR__
{
    public:
        sR2__(YCS name, YCI t, YCQR a1, YCQR a2) : sR__(name, t, a1)
        {   
            pars_.push_back(a2);
        }

        YSG copy_gate() const { return std::make_shared<sR2__>(*this); };
};


/** Ry(angle_ry) * Rz(angle_rz) */
class Rc__ : public sR2__
{
    public:
        Rc__(YCI t, YCQR angle_rz, YCQR angle_ry) : sR2__(name_shared_, t, angle_rz, angle_ry)
        { 
            u2_ = YGV::mRc(angle_rz, angle_ry); 
            tex_name_ = "R_c";
        }
        YSG copy_gate() const { return std::make_shared<Rc__>(*this); };

    public:
        const static std::string name_shared_;
};


/**
 * Rc(angle_rz, angle_ry) that takes action only if is put inside a gadget;
*/
class Rc_gadget__ : public SQGate__
{
    public:
        Rc_gadget__(YCI t, YCVQ angles_y, YCVQ angles_z, YCVI r_counter) : SQGate__(name_shared_, t)
        { 
            cs_counter_ = r_counter;
            angles_y_ = angles_y;
            angles_z_ = angles_z;

            // check the number of angles:
            if(angles_z_.size() != angles_y_.size())
            {
                throw std::string("Error in Rc_gadget__: the number of angles are not equal.");
            }
            flag_activated_ = false;
            tex_name_ = "R_c";
        }
        Rc_gadget__(const Rc_gadget__& oo) : SQGate__(oo)
        {
            cs_counter_      = oo.cs_counter_;
            angles_y_       = oo.angles_y_;
            angles_z_       = oo.angles_z_;
            flag_activated_ = oo.flag_activated_;
            id_in_counter_  = oo.id_in_counter_;
        }
        void activate_gadget(YCI id_in_counter, YCI N_mult)
        {
            flag_activated_ = true;
            id_in_counter_ = id_in_counter;
            if((N_mult-id_in_counter_) >= angles_z_.size())
            {
                std::cout << "In Rc_gadget: index = " << N_mult-id_in_counter_ << std::endl;
                throw std::string("Error in Rc_gadget__: the index is too large");
            }

            // std::cout << id_in_counter_ << std::endl;
            // std::cout << "activation, ay, az: " << 
            //     angles_y_[N_mult-id_in_counter_] << ", " << angles_z_[N_mult-id_in_counter_] << std::endl;

            u2_ = YGV::mRc(angles_z_[N_mult-id_in_counter_], angles_y_[N_mult-id_in_counter_]);
        }
        void deactivate_gadget()
        {
            flag_activated_ = false;
        }
        void generate(Qureg& oc)
        {
            if(flag_activated_)
            {
                int nc = cs_counter_.size();
                YVIv cs_unit_total, cs_zero_total;

                std::vector<short> binArray(nc);
                YMATH::intToBinary(id_in_counter_, binArray);
                for(unsigned id_bit = 0; id_bit < nc; id_bit++)
                    if(binArray[nc - id_bit - 1] == 1)
                        cs_unit_total.push_back(cs_counter_[id_bit]);
                    else 
                        cs_zero_total.push_back(cs_counter_[id_bit]);
                cs_unit_total.insert(cs_unit_total.end(), cs_unit_.begin(), cs_unit_.end());
                cs_zero_total.insert(cs_zero_total.end(), cs_zero_.begin(), cs_zero_.end());

                // std::cout << "\ngeneration: " << id_in_counter_ << std::endl;
                // YMIX::print(cs_unit_total);
                // YMIX::print(cs_zero_total);

                mc_st_u(oc, ts_[0], cs_unit_total, cs_zero_total, u2_);
            }
        }
        YSG copy_gate() const { return std::make_shared<Rc_gadget__>(*this); };

        void write_tex(
            std::vector<std::vector<std::string>>& tex_lines, 
            const int64_t& id_layer,
            YCU nq
        ){
            if(flag_activated_)
            {
                SQGate__::write_tex(tex_lines, id_layer, nq);
            }
        }

    public:
        // qubits by which the gate will be controlled in compression gadget.
        YVIv cs_counter_; 

        // integer index indicating which angle_y and angle_z to use in the gate.
        int id_in_counter_;

        /**
         * \p angles_y_[id_in_counter_] and \p angles_z_[id_in_counter_] are used in the gate 
         * when the register \p cs_counter_ encodes the integer \p id_in_counter_.
        */
        YVQv angles_y_, angles_z_; 

        /**
         * If true, then gate is applied.
        */
        bool flag_activated_;

        const static std::string name_shared_;
};





#endif