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
    unsigned nFor;      // number of foraging events for i
    unsigned nS;        // number of scroungers at current event
    flt p;              // prob for i to use action 1
    flt R;              // perceived reward by i
    flt delt;           // TD error for i
    flt Q0;             // estimated value of action 0 for i
    flt Q1;             // estimated value of action 1 for i
};

//************************ class ActValGroup *****************************

// This class sets up and simulates the action-value (Sarsa) learning method
// for a producer-scrounger game.

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
    using rand_discr = std::discrete_distribution<int>;
    ActValGroup(unsigned a_N,
                unsigned a_KR,
                unsigned a_nScr,
                unsigned a_T,
                unsigned a_T0hist,
                flt a_prwr,
                flt a_V0,
                flt a_V1,
                flt a_V2,
                flt a_prs,
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
    unsigned nScr;   // maximum number of scroungers at foraging event
    unsigned T;      // total number of rounds for group
    unsigned T0hist; // round to start collecting learning history
    flt prwr;        // probability to "rewire edge" in graph of connections
    flt V0;          // payoff parameter
    flt V1;          // payoff parameter
    flt V2;          // payoff parameter
    flt prs;         // factor for probability of success for producer
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
    unsigned a_nScr,
    unsigned a_T,
    unsigned a_T0hist,
    flt a_prwr,
    flt a_V0,
    flt a_V1,
    flt a_V2,
    flt a_prs,
    flt a_pmax,
    flt a_bet,
    const vph_type& a_memb,
    rand_eng& a_eng,
    bool a_lhist) :
    N{a_N},
    KR{a_KR},
    nScr{a_nScr},
    T{a_T},
    T0hist{a_T0hist},
    prwr{a_prwr},
    V0{a_V0},
    V1{a_V1},
    V2{a_V2},
    prs{a_prs},
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
    rand_uni Vm(0.5, 1.5);
    // set payoff values to V0 at start of generation
    for (auto& m : memb) {
        m.payoff = V0;
    }
    int nscr = nScr;
    int Kr = KR;
    // unsigned K = 2*KR;
    for (unsigned tstep = 0; tstep < T; ++tstep) {
        // Producing prob and indicator
        v_type pP(N, 0);
        ui_type u(N, 0);
        // number of producers among connected individuals
        i_type kP(N, 0);
        // rewards from current time step
        v_type R(N, 0);
        // in each time step, each group member decides whether to Produce or
        // to Scrounge
        for (unsigned i = 0; i < N; ++i) {
            pP[i] = memb[i].p1(bet);
            u[i] = (uni(eng) < Clamp(pP[i])) ? 1:0;
            // u[i] = (uni(eng) < pP[i]) ? 1:0;
            memb[i].u = u[i];
            memb[i].nS = 0;
        }
        // compute polarisation
        Fv[tstep] = F();
        // work out kP (number of producers connected to)
        for (int i = 0; i < N; ++i) {
            // // NOTE: for code without boost graph
            // for (int jd = -Kr; jd <= Kr; ++jd) {
            //     if (jd != 0) {
            //         int j = (i + jd + N) % N;
            //         if (u[j] == 1) ++kP[i];
            //     }
            // }
            if (g.complete) {
                // all group members are connected
                for (int j = 0; j < N; ++j) {
                    if (j != i) {
                        if (u[j] == 1) ++kP[i];
                    }
                }
            } else {
                adj_iter_range itr = adjacent_vertices(i, g.gr);
                for (auto cur = itr.first; cur != itr.second; ++cur) {
                    int j = *cur;
                    if (u[j] == 1) ++kP[i];
                }
            }
            memb[i].kP = kP[i];
        }
        // run through individuals and process foraging events (each event where
        // a producer finds food is treated as separate in time, in the sense
        // that scroungers can in principle be present at separate events during
        // a time step)
        for (unsigned i = 0; i < N; ++i) {
            // check if i is producer and found food
            if (u[i] == 1 && uni(eng) < prs*memb[i].q) {
                flt vm = Vm(eng);
                ++memb[i].nFor;
                // if (kP[i] < K) { // there are potential scroungers
                if (kP[i] < g.k[i]) { // there are potential scroungers
                    // actual number of scroungers
                    // int nS = std::min(K - kP[i], nScr);
                    int nS = std::min(g.k[i] - kP[i], nscr);
                    // i_type iPot(K - kP[i], 0);
                    i_type iPot(g.k[i] - kP[i], 0);
                    unsigned cnt = 0;
                    // for (int jd = -Kr; jd <= Kr; ++jd) {
                    //     if (jd != 0) {
                    //         int j = (i + jd + N) % N;
                    //         if (u[j] == 0) { // a scrounger
                    //             iPot[cnt] = j;
                    //             ++cnt;
                    //         }
                    //     }
                    // }
                    if (g.complete) {
                        // all group members are connected
                        for (int j = 0; j < N; ++j) {
                            if (u[j] == 0 && j != i) { // a scrounger
                                iPot[cnt] = j;
                                ++cnt;
                            }
                        }
                    } else {
                        adj_iter_range itr = adjacent_vertices(i, g.gr);
                        for (auto cur = itr.first; cur != itr.second; ++cur) {
                            int j = *cur;
                            if (u[j] == 0) { // a scrounger
                                iPot[cnt] = j;
                                ++cnt;
                            }
                        }
                    }
                    // randomly shuffle indices of potential scroungers
                    std::shuffle(iPot.begin(), iPot.end(), eng);
                    // sample nS actual scroungers from the potential ones
                    i_type iAct;
                    iAct.reserve(nS);
                    std::sample(iPot.begin(), iPot.end(),
                        std::back_inserter(iAct), nS, eng);
                    // update rewards and payoffs
                    // first the producer
                    R[i] += vm*(V1 + V2/(1 + nS));
                    memb[i].payoff += vm*(V1 + V2/(1 + nS));
                    memb[i].nS = nS;
                    // then the nS scroungers
                    for (unsigned l = 0; l < nS; ++l) {
                        int j = iAct[l];
                        ++memb[j].nFor;
                        memb[j].nS = nS;
                        R[j] += memb[j].v*vm*(V2/(1 + nS));
                        memb[j].payoff += vm*(V2/(1 + nS));
                    }
                } else { // no connected scroungers; only producer present
                    memb[i].nS = 0;
                    R[i] += vm*(V1 + V2);
                    memb[i].payoff += vm*(V1 + V2);
                }
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
    flt varp = ssp/N - av_p*av_p;
    return varp/(av_p*(1 - av_p));
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
    st.nFor = mi.nFor;
    st.nS = mi.nS;
    st.p = mi.p;
    st.R = mi.R;
    st.delt = delt;
    st.Q0 = mi.Q0;
    st.Q1 = mi.Q1;
    stat.push_back(st);
}

#endif // AVGROUP_HPP
