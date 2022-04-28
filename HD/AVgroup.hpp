#ifndef AVGROUP_HPP
#define AVGROUP_HPP

#include "SWgraph.hpp"
#include <vector>
#include <random>
#include <cmath>
#include <algorithm>

// The Speci program runs simulations of learning in a social network
// Copyright (C) 2022  Olof Leimar
// See Readme.md for copyright notice

//************************* struct LearnStat ******************************

// This struct stores data on an interaction between two group members (i and
// j); a sequence of structs can be saved as a record of the interaction
// history in the group

template<typename PhenType>
struct LearnStat {
    using phen_type = PhenType;
    using flt = typename phen_type::flt;
    unsigned gnum;      // group number
    unsigned tstep;     // time step of interaction in group
    unsigned i;         // inum for individual
    unsigned u;         // action (0 or 1)
    flt p;              // prob for i to use action 1
    flt R;              // perceived reward by i
    flt delt;           // TD error for i
    flt Q0;             // estimated value of action 0 for i
    flt Q1;             // estimated value of action 1 for i
};

//************************ class ActValGroup *****************************

// This class sets up and simulates the Sarsa learning method for a
// Hawk-Dove game.

// The class deals with the interactions in one group, over the time steps
// during one generation.

template<typename PhenType>
class ActValGroup {
public:
    using phen_type = PhenType;
    using vph_type = std::vector<phen_type>;
    using flt = typename phen_type::flt;
    using v_type = std::vector<flt>;
    using ui_type = std::vector<unsigned>;
    using stat_type = LearnStat<phen_type>;
    using vs_type = std::vector<stat_type>;
    using gcont_type = SWGraph<flt>;
    using graph_type = typename gcont_type::graph_type;
    using adj_iter_type = typename gcont_type::adj_iter_type;
    using adj_iter_range = typename gcont_type::adj_iter_range;
    using rand_eng = std::mt19937;
    using rand_uni = std::uniform_real_distribution<flt>;
    using rand_int = std::uniform_int_distribution<int>;
    using rand_norm = std::normal_distribution<flt>;
    ActValGroup(unsigned a_N,
                unsigned a_KR,
                unsigned a_T,
                unsigned a_T0hist,
                flt a_prwr,
                flt a_V0,
                flt a_V,
                flt a_C,
                flt a_pmax,
                flt a_bet,
                const vph_type& a_memb,
                rand_eng& a_eng,
                bool a_lhist = false);
    const vph_type& Get_memb() const { return memb; }
    const vs_type& Get_stat() const { return stat; }
    const v_type& Get_F() const { return Fv; }
    void Interact();

private:
    // prob i wins over j in Hawk-Hawk fight
    flt phi(flt qi, flt qj) { return (1.0 + qi - qj)/2; }
    flt Clamp(flt p)
        {return (p > pmax) ? pmax : ((p < 1 - pmax) ? 1 - pmax : p); }
    void Add_lstat(unsigned tstep, unsigned i, flt delt);
    flt F();
    unsigned N;      // group size
    unsigned KR;     // the degree of the graph is 2*KR
    unsigned T;      // total number of rounds for group
    unsigned T0hist; // round to start collecting learning history
    flt prwr;        // probability to "rewire edge" in graph of connections
    flt V0;          // payoff parameter
    flt V;           // payoff parameter
    flt C;           // payoff parameter
    flt pmax;        // maximum value for probability to use action 1
    flt bet;         // parameter for soft-max function
    vph_type memb;   // phenotypes of members of the group
    rand_eng& eng;   // reference to random number engine
    gcont_type g;    // contains graph of connections between group members
    v_type Fv;       // storage for polarisation indices
    bool lhist;      // whether to collect learning history
    vs_type stat;    // learning statistics
};

template<typename PhenType>
ActValGroup<PhenType>::ActValGroup(unsigned a_N,
    unsigned a_KR,
    unsigned a_T,
    unsigned a_T0hist,
    flt a_prwr,
    flt a_V0,
    flt a_V,
    flt a_C,
    flt a_pmax,
    flt a_bet,
    const vph_type& a_memb,
    rand_eng& a_eng,
    bool a_lhist) :
    N{a_N},
    KR{a_KR},
    T{a_T},
    T0hist{a_T0hist},
    prwr{a_prwr},
    V0{a_V0},
    V{a_V},
    C{a_C},
    pmax{a_pmax},
    bet{a_bet},
    memb{a_memb},
    eng{a_eng},
    g(N, KR, prwr, eng),
    lhist{a_lhist}
{
    // introduce initial random variation in Q0 and Q1
    rand_uni dQ(-0.01, 0.01);
    for (unsigned i = 0; i < N; ++i) {
        memb[i].Q0 += dQ(eng);
        memb[i].Q1 += dQ(eng);
    }
    Fv.reserve(T);
    if (lhist) {
        stat.reserve(T);
    }
}

template<typename PhenType>
void ActValGroup<PhenType>::Interact()
{
    // // NOTE: temporary code to write graph to file
    // g.WriteGraphML("Run_g21_deg04_prw01.graphml");
    rand_uni uni(0, 1);
    rand_int uri(0, N - 1);
    // NOTE: for code without boost graph
    // rand_int rci(-KR, KR);
    rand_int incr(0, 2*KR - 1); // default increment to out edge iterator
    // set payoff values to zero at start of generation
    for (auto& m : memb) {
        m.payoff = 0;
    }
    unsigned Th = T/2; // each HD contest involves 2 individuals
    for (unsigned tstep2 = 0; tstep2 < N*Th; ++tstep2) {
        // let tstep be equal to expected number of contests per individual
        unsigned tstep = 2*(tstep2/N);
        // select random individual
        int i = uri(eng);
        int j = 0;
        int l = 0;
        // select random connected opponent
        // // NOTE: for code without boost graph
        // while (l == 0) {
        //     l = rci(eng);
        // }
        // j = (i + l + N) % N;
        if (g.complete) {
            // all individuals in group are connected
            j = i;
            while (j == i) {
                j = uri(eng);
            }
        } else {
            if (g.k[i] == 2*KR) {
                // use default increment for out edge iterator
                l = incr(eng);
            } else if (g.k[i] > 0) {
                // construct appropriate random increment
                rand_int incr2(0, g.k[i] - 1);
                l = incr2(eng);
            } else {
                // no connected opponent; move on to next i
                continue;
            }
            adj_iter_range it_range = adjacent_vertices(i, g.gr);
            adj_iter_type fst = it_range.first;
            adj_iter_type opp = fst + l;
            j = *opp; // this is the random neighbour opponent
        }
        // go ahead with contest
        phen_type& mi = memb[i];
        phen_type& mj = memb[j];
        // probability using Hawk
        flt pij = mi.p1(bet);
        flt pji = mj.p1(bet);
        // compute polarisation every second contest per individual
        if (tstep2 % N == 0) {
            Fv.push_back(F());
        }
        unsigned ui = (uni(eng) < Clamp(pij)) ? 1:0;
        unsigned uj = (uni(eng) < Clamp(pji)) ? 1:0;
        // unsigned ui = (uni(eng) < pij) ? 1:0;
        // unsigned uj = (uni(eng) < pji) ? 1:0;
        mi.u = ui;
        mj.u = uj;
        mi.nInts += 1;
        mj.nInts += 1;
        bool iwin = (ui == 1);
        if (ui*uj > 0) { // Hawk-Hawk fight
            mi.nHH += 1;
            mj.nHH += 1;
            iwin = (uni(eng) < phi(mi.q, mj.q));
        } else if (ui == 0 && uj == 0) { // Dove-Dove
            iwin = (uni(eng) < 0.5);
        }
        flt Rij = 0;
        flt Rji = 0;
        if (ui == 1) {
            if (uj == 1) { // Hawk-Hawk
                if (iwin) {
                    Rij = mi.v*V;
                    Rji = -C;
                    // Rij = (mi.v*V - C)/2;
                    // Rji = (mj.v*V - C)/2;
                } else {
                    Rij = -C;
                    Rji = mj.v*V;
                    // Rij = (mi.v*V - C)/2;
                    // Rji = (mj.v*V - C)/2;
                }
            } else {
                Rij = mi.v*V;
                Rji = 0;
            }
        } else {
            if (uj == 1) {
                Rij = 0;
                Rji = mj.v*V;
            } else { // Dove-Dove
                if (iwin) {
                    Rij = mi.v*V;
                    Rji = 0;
                    // Rij = mi.v*V/2;
                    // Rji = mj.v*V/2;
                } else {
                    Rij = 0;
                    Rji = mj.v*V;
                    // Rij = mi.v*V/2;
                    // Rji = mj.v*V/2;
                }
            }
        }
        // perform payoff increments
        mi.payoff += V0;
        mj.payoff += V0;
        if (ui == 1) {
            if (uj == 1) {
                mi.payoff += iwin ? V : -C;
                mj.payoff += iwin ? -C : V;
            } else {
                mi.payoff += V;
            }
        } else {
            if (uj == 1) {
                mj.payoff += V;
            } else {
                mi.payoff += iwin ? V : 0;
                mj.payoff += iwin ? 0 : V;
            }
        }
        // update learning parameters
        flt delt = 0;
        if (ui == 0) {
            delt = Rij - mi.Q0;
            mi.Q0 += mi.alph*delt;
        } else {
            delt = Rij - mi.Q1;
            mi.Q1 += mi.alph*delt;
        }
        mi.R = Rij;
        if (lhist && tstep > T0hist) {
            Add_lstat(tstep, i, delt);
        }
        // update learning parameters
        if (uj == 0) {
            delt = Rji - mj.Q0;
            mj.Q0 += mj.alph*delt;
        } else {
            delt = Rji - mj.Q1;
            mj.Q1 += mj.alph*delt;
        }
        mj.R = Rji;
        if (lhist && tstep > T0hist) {
            Add_lstat(tstep, j, delt);
        }
    }
    // scale payoff to be per interaction
    for (auto& m : memb) {
        if (m.nInts > 0) {
            m.payoff /= m.nInts;
        }
    }
}

template<typename PhenType>
typename ActValGroup<PhenType>::flt ActValGroup<PhenType>::F()
{
    flt av_p = 0;
    for (auto& m : memb) {
        av_p += m.p;
    }
    av_p /= N;
    flt ssp = 0;
    for (auto& m : memb) {
        ssp += m.p*m.p;
    }
    ssp /= N;
    flt varp = ssp - av_p*av_p;
    flt F = varp/(av_p*(1 - av_p));
    return F;
}

template<typename PhenType>
void ActValGroup<PhenType>::Add_lstat(unsigned tstep, unsigned i, flt delt)
{
    phen_type& mi = memb[i];
    stat_type st;
    st.gnum = mi.gnum;          // group number
    st.tstep = tstep;           // time step (round) for group
    st.i = mi.inum;
    st.u = mi.u;
    st.p = mi.p;
    st.R = mi.R;
    st.delt = delt;
    st.Q0 = mi.Q0;
    st.Q1 = mi.Q1;
    stat.push_back(st);
}

#endif // AVGROUP_HPP
