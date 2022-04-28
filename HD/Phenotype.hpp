#ifndef PHENOTYPE_HPP
#define PHENOTYPE_HPP

#include <string>
#include <array>
#include <ostream>
#include <istream>
#include <cmath>

// The Speci program runs simulations of learning in a social network
// Copyright (C) 2022  Olof Leimar
// See Readme.md for copyright notice

//************************* struct Phenotype ******************************

// Assumptions about GenType:
// types:
//   val_type
// member functions:
//   val_type Value()

// In addition to the genotypic trait values, consisting of Qini, alph, v, this
// class also stores the individual quality q, the average per round payoff,
// the current value of the learning parameters (when saved this will be the
// value after the specified number of rounds during a generation), individual
// number in group, the group number, the subpopulation number, the
// individual's sex, and whether it is alive.

template<typename GenType>
struct Phenotype {
// public:
    using flt = double;
    using v_type = std::vector<flt>;
    using gen_type = GenType;
    using val_type = typename gen_type::val_type;
    Phenotype(flt a_Qini,
        flt a_alph,
        flt a_v,
        flt a_q,
        flt a_payoff,
        unsigned a_inum,
        unsigned a_gnum,
        bool a_female,
        bool a_alive) :
        Qini{a_Qini},
        alph{a_alph},
        v{a_v},
        q{a_q},
        Q0{Qini},
        Q1{Qini},
        p{0},
        payoff{a_payoff},
        R{0},
        u{0},
        inum{a_inum},
        gnum{a_gnum},
        female{a_female},
        alive{a_alive} {}
    Phenotype(const gen_type& gt = gen_type()) { Assign(gt); }
    void Assign(const gen_type& gt);
    void Set_inum(unsigned a_inum);
    // probability of choosing Hawk
    flt p1(flt bet) { return p = 1/(1 + std::exp(-bet*(Q1 - Q0))); }
    bool Female() const { return female; }
    // public data members
    flt Qini;    // parameter for estimated value at start of a generation
    flt alph;    // learning rate
    flt v;       // perceived value factor
    flt q;       // individual quality (fighting ability)
    flt Q0;      // estimated value of action u = 0 (Dove)
    flt Q1;      // estimated value of action u = 1 (Hawk)
    flt p;       // probability to play Hawk
    flt payoff;  // the (fitness) payoff per round
    flt R;       // reward in recent time step
    unsigned u;      // action used (0 is Dove, 1 is Hawk)
    unsigned nInts;  // number of Hawk-Dove interactions
    unsigned nHH;    // number of Hawk-Hawk interactions
    unsigned inum;   // individual number
    unsigned gnum;   // group number
    bool female;
    bool alive;
};

template<typename GenType>
void Phenotype<GenType>::Assign(const gen_type& gt)
{
    val_type val = gt.Value();
    // assume val is a vector with components corresponding to the traits Qini,
    // alph, and v
    Qini = val[0];
    alph = val[1];
    v = val[2];
    q = 0;
    Q0 = Qini;
    Q1 = Qini;
    p = 0.5;
    payoff = 0;
    R = 0;
    u = 0;
    nInts = 0;
    nHH = 0;
    inum = 0;
    gnum = 0;
    female = true;
    alive = true;
}

template<typename GenType>
void Phenotype<GenType>::Set_inum(unsigned a_inum)
{
    inum = a_inum;
}

#endif // PHENOTYPE_HPP
