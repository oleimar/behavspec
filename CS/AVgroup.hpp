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

// This class sets up and simulates action-value learning (Sarsa) for a
// caller-satellite game.

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
    using i_type = std::vector<int>;
    using stat_type = LearnStat<phen_type>;
    using vs_type = std::vector<stat_type>;
    using gcont_type = SWGraph<flt>;
    using graph_type = typename gcont_type::graph_type;
    using adj_iter_type = typename gcont_type::adj_iter_type;
    using adj_iter_range = typename gcont_type::adj_iter_range;
    using rand_eng = std::mt19937;
    using rand_uni = std::uniform_real_distribution<flt>;
    using rand_int = std::uniform_int_distribution<int>;
    using rand_poiss = std::poisson_distribution<int>;
    using rand_discr = std::discrete_distribution<int>;
    ActValGroup(unsigned a_N,
                unsigned a_KR,
                unsigned a_T,
                unsigned a_T0hist,
                flt a_prwr,
                flt a_V0,
                flt a_V1,
                flt a_gam0,
                flt a_f,
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
    flt Clamp(flt p)
        {return (p > pmax) ? pmax : ((p < 1 - pmax) ? 1 - pmax : p); }
    flt F();
    void Add_lstat(unsigned tstep, unsigned i, flt delt);
    unsigned N;      // group size
    unsigned KR;     // the degree of the graph is 2*KR
    unsigned T;      // total number of rounds for group
    unsigned T0hist; // round to start collecting learning history
    flt prwr;        // probability to "rewire edge" in graph of connections
    flt V0;          // payoff parameter
    flt V1;          // payoff parameter
    flt gam0;        // payoff parameter
    flt f;           // female arrival rate
    flt pmax;        // maximum value for probability to use action
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
    flt a_V1,
    flt a_gam0,
    flt a_f,
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
    V1{a_V1},
    gam0{a_gam0},
    f{a_f},
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
    Fv.resize(T, 0);
    if (lhist) {
        stat.reserve(T);
    }
}

template<typename PhenType>
void ActValGroup<PhenType>::Interact()
{
    rand_uni uni(0, 1);
    rand_int uri(0, N - 1);
    // set payoff values to V0 at start of generation
    for (auto& m : memb) {
        m.payoff = V0;
    }
    // NOTE: call_pen is used to indicate how rewards are implemented
    bool call_pen = false;
    for (unsigned tstep = 0; tstep < T; ++tstep) {
        // calling prob and indicator
        v_type pC(N, 0);
        ui_type u(N, 0);
        // number of callers among connected males
        i_type kC(N, 0);
        // call strengths
        v_type s(N, 0);
        // perceived neighbour aggression (only used if call_pen is true)
        v_type pna(N, 0);
        // rewards from current round
        v_type R(N, 0);
        // in each time step, each group member decides whether to Call or
        // to act as Satellite
        unsigned nC = 0;
        for (unsigned i = 0; i < N; ++i) {
            pC[i] = memb[i].p1(bet);
            u[i] = (uni(eng) < Clamp(pC[i])) ? 1:0;
            // u[i] = (uni(eng) < pC[i]) ? 1:0;
            memb[i].u = u[i];
            nC += u[i];
        }
        // compute polarisation
        Fv[tstep] =F();
        // work out connected callers, call strength, and pna[i]
        flt sums = 0; // sum of all call strengths
        for (int i = 0; i < N; ++i) {
            // get kC both for callers and satellites
            if (g.complete) {
                // all group members are connected
                for (int j = 0; j < N; ++j) {
                    if (j != i) {
                        if (u[j] == 1) ++kC[i];
                    }
                }
            } else {
                adj_iter_range itr = adjacent_vertices(i, g.gr);
                for (auto cur = itr.first; cur != itr.second; ++cur) {
                    int j = *cur;
                    if (u[j] == 1) ++kC[i];
                }
            }
            memb[i].kC = kC[i];
            if (u[i] == 1) { // set call strength for caller
                s[i] = 1 - (gam0*kC[i])/g.k[i];
                sums += s[i];
            }
            memb[i].s = s[i];
            // perceived neighbour aggression
            if (u[i] == 1) {
                pna[i] = (memb[i].ppf*kC[i])/g.k[i];
            }
            memb[i].pna = pna[i];
        }
        if (call_pen) {
            for (int i = 0; i < N; ++i) {
                if (u[i] == 1) {
                    R[i] = 1 - pna[i];
                }
            }
        }
        // number of females arriving this time step
        rand_poiss poi(f*sums);
        unsigned nf = poi(eng);
        for (unsigned ii = 0; ii < nf; ++ii) {
            rand_discr dscr(s.begin(), s.end());
            int i = dscr(eng);
            int kSi = g.k[i] - kC[i];
            if (kSi > 0) { // there are satellites
                i_type ind(kSi + 1, 0);
                v_type wei(kSi + 1, 0);
                unsigned cnt = 0;
                ind[cnt] = i;   // the caller
                wei[cnt] = 1.0; // caller has weight 1
                if (g.complete) {
                    // all group members are connected
                    for (int j = 0; j < N; ++j) {
                        if (u[j] == 0 && j != i) { // a satellite
                            ++cnt;
                            ind[cnt] = j;
                            wei[cnt] = 1.0/kSi;
                        }
                    }
                } else {
                    adj_iter_range itr = adjacent_vertices(i, g.gr);
                    for (auto cur = itr.first; cur != itr.second; ++cur) {
                        int j = *cur;
                        if (u[j] == 0) { // a satellite
                            ++cnt;
                            ind[cnt] = j;
                            wei[cnt] = 1.0/kSi;
                        }
                    }
                }
                rand_discr dml(wei.begin(), wei.end());
                int jn = dml(eng);
                int j = ind[jn]; // this is the male that mates
                if (!call_pen) {
                    R[j] += V1;
                }
                memb[j].payoff += V1;
            } else { // female mates with calling male
                if (!call_pen) {
                    R[i] += V1;
                }
                memb[i].payoff += V1;
            }
        }
        // update learning parameters
        for (unsigned i = 0; i < N; ++i) {
            flt delt = 0;
            phen_type& mi = memb[i];
            if (u[i] == 0) {
                delt = R[i] - mi.Q0;
                mi.Q0 += mi.alph*delt;
            } else {
                delt = R[i] - mi.Q1;
                mi.Q1 += mi.alph*delt;
            }
            mi.R = R[i];
            if (lhist && tstep > T0hist) {
                Add_lstat(tstep, i, delt);
            }
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
