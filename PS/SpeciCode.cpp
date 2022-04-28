#include "cpptoml.h"   // to read input parameters from TOML file
#include "SpeciCode.hpp"
#include "hdf5code.hpp"
#include "Utils.hpp"
#include <algorithm>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>

// The Speci program runs simulations of learning in a social network
// Copyright (C) 2022  Olof Leimar
// See Readme.md for copyright notice

#ifdef PARA_RUN
#include <omp.h>
#endif

//************************** Read and ReadArr ****************************

// convenience functions to read from TOML input file

// this template function can be used for any type of single value
template<typename T>
void Get(std::shared_ptr<cpptoml::table> infile,
         T& value, const std::string& name)
{
    auto val = infile->get_as<T>(name);
    if (val) {
        value = *val;
    } else {
        std::cerr << "Read failed for identifier " << name << "\n";
    }
}

// this template function can be used for a vector or array (but there is no
// checking how many elements are read)
template<typename It>
void GetArr(std::shared_ptr<cpptoml::table> infile,
            It beg, const std::string& name)
{
    using valtype = typename std::iterator_traits<It>::value_type;
    auto vp = infile->get_array_of<valtype>(name);
    if (vp) {
        std::copy(vp->begin(), vp->end(), beg);
    } else {
        std::cerr << "Read failed for identifier " << name << "\n";
    }
}


//************************** class SpeciInpData ****************************

SpeciInpData::SpeciInpData(const char* filename) :
      OK(false)
{
    auto idat = cpptoml::parse_file(filename);
    Get(idat, max_num_thrds, "max_num_thrds");
    Get(idat, num_loci, "num_loci");
    Get(idat, ng, "ng");
    Get(idat, N, "N");
    Get(idat, KR, "KR");
    Get(idat, nScr, "nScr");
    Get(idat, T, "T");
    Get(idat, T0hist, "T0hist");
    Get(idat, numgen, "numgen");
    Get(idat, prwr, "prwr");
    Get(idat, V0, "V0");
    Get(idat, V1, "V1");
    Get(idat, V2n, "V2n");
    Get(idat, V2x, "V2x");
    Get(idat, qn, "qn");
    Get(idat, qx, "qx");
    Get(idat, prs, "prs");
    Get(idat, pmax, "pmax");
    Get(idat, bet, "bet");
    mut_rate.resize(num_loci);
    GetArr(idat, mut_rate.begin(), "mut_rate");
    SD.resize(num_loci);
    GetArr(idat, SD.begin(), "SD");
    max_val.resize(num_loci);
    GetArr(idat, max_val.begin(), "max_val");
    min_val.resize(num_loci);
    GetArr(idat, min_val.begin(), "min_val");
    rho.resize(num_loci);
    GetArr(idat, rho.begin(), "rho");
    Get(idat, learn_hist, "learn_hist");
    Get(idat, read_from_file, "read_from_file");
    if (read_from_file) {
        Get(idat, h5InName, "h5InName");
    } else {
        all0.resize(num_loci);
        GetArr(idat, all0.begin(), "all0");
    }
    Get(idat, h5OutName, "h5OutName");
    if (learn_hist) {
        Get(idat, h5HistName, "h5HistName");
    }
    InpName = std::string(filename);
    OK = true;
}


//****************************** Class Speci *****************************

Speci::Speci(const SpeciInpData& sid) :
    id{sid},
    num_loci{id.num_loci},
    ng{id.ng},
    N{id.N},
    KR{id.KR},
    nScr{id.nScr},
    Ntot{ng*N},
    T{id.T},
    T0hist{id.T0hist},
    numgen{id.numgen},
    prwr{static_cast<flt>(id.prwr)},
    V0{static_cast<flt>(id.V0)},
    V1{static_cast<flt>(id.V1)},
    V2n{static_cast<flt>(id.V2n)},
    V2x{static_cast<flt>(id.V2x)},
    qn{static_cast<flt>(id.qn)},
    qx{static_cast<flt>(id.qx)},
    prs{static_cast<flt>(id.prs)},
    pmax{static_cast<flt>(id.pmax)},
    bet{static_cast<flt>(id.bet)},
    learn_hist{id.learn_hist},
    num_thrds{1}
{
    // decide on number of threads for parallel processing
#ifdef PARA_RUN
    num_thrds = omp_get_max_threads();
    if (num_thrds > id.max_num_thrds) num_thrds = id.max_num_thrds;
    std::cout << "Number of threads: " << num_thrds << '\n';
#endif
    // generate one seed for each thread
    sds.resize(num_thrds);
    std::random_device rd;
    for (unsigned i = 0; i < num_thrds; ++i) {
        sds[i] = rd();
        vre.push_back(rand_eng(sds[i]));
    }

    // Note concerning thread safety: in order to avoid possible problems with
    // multiple threads, the std::vector container pop is allocated once and
    // for all here, and thread-local data are then copied into position in
    // pop (thus avoiding potentially unsafe push_back and insert).

    // create Ntot "placeholder individuals" in population
    gam_type gam(num_loci);
    ind_type indi(gam);
    pop.resize(Ntot, indi);
    // learning history stats
    if (learn_hist) {
        stat.reserve(ng*T);
    }

    // check if population data should be read from file
    if (id.read_from_file) {
        // Read_pop(id.InName);
        h5_read_pop(id.h5InName);
    } else {
        // construct all individuals as essentially the same
        gam_type gam(num_loci); // starting gamete
        for (unsigned l = 0; l < num_loci; ++l) {
            gam.gamdat[l] = static_cast<flt>(id.all0[l]);
        }
        unsigned j = 0;
        for (unsigned gn = 0; gn < ng; ++gn) { // groups
            for (unsigned i = 0; i < N; ++i) { // inds in group
                ind_type ind(gam);
                ind.phenotype.gnum = gn; // set group number
                ind.phenotype.Set_inum(i); // set individual number
                pop[j++] = ind;
            }
        }
    }
    avF.resize(T, 0);
}

void Speci::Run()
{
    Timer timer(std::cout);
    timer.Start();
    ProgressBar PrBar(std::cout, numgen);
    // set up  "global" mutation record, with engine and parameters controlling
    // mutation, segregation and recombination
    mut_rec_type mr0(vre[0], num_loci);
    for (unsigned l = 0; l < num_loci; ++l) {
        mr0.mut_rate[l] = static_cast<flt>(id.mut_rate[l]);
        mr0.SD[l] = static_cast<flt>(id.SD[l]);
        mr0.max_val[l] = static_cast<flt>(id.max_val[l]);
        mr0.min_val[l] = static_cast<flt>(id.min_val[l]);
        mr0.rho[l] = static_cast<flt>(id.rho[l]);
    }
    // run through generations
    for (unsigned gen = 0; gen < numgen; ++gen) {
        // use parallel for processing over the interactions in groups
#pragma omp parallel for num_threads(num_thrds)
        for (unsigned gn = 0; gn < ng; ++gn) {
#ifdef PARA_RUN
            int threadn = omp_get_thread_num();
#else
            int threadn = 0;
#endif
            // thread-local random number engine
            rand_eng& eng = vre[threadn];
            rand_uni uni(0, 1);
            rand_uni unq(qn, qx);
            rand_uni unV2(V2n, V2x);

            // thread-local container for group phenotypes
            vph_type gph;
            gph.reserve(N);
            for (unsigned i = 0; i < N; ++i) {
                // copy phenotypes from population to group
                gph.push_back(pop[gn*N + i].phenotype);
            }

            // assign (random) quality q to group members
            for (unsigned i = 0; i < N; ++i) {
                gph[i].q = unq(eng);
            }
            // assign (random) value of V2
            V2 = unV2(eng);

            // get history only for single threaded and final generation
            bool lhist = false;
            if (num_thrds == 1 && gen == numgen - 1 && learn_hist) {
                lhist = true;
            }

            // set up group interaction
            avg_type avg(N, KR, nScr, T, T0hist, prwr, V0, V1, V2, prs,
                pmax, bet, gph, eng, lhist);
            avg.Interact();
            // get resulting phenotype of group members
            const vph_type& memb = avg.Get_memb();
            // copy individuals from group back into (global) container
            for (unsigned i = 0; i < N; ++i) {
                pop[gn*N + i].phenotype = memb[i];
            }
            // update global avF
            const v_type& gF = avg.Get_F();
            for (unsigned tp = 0; tp < T; ++tp) {
                avF[tp] += gF[tp];
            }
            if (lhist) {
                // append history from the group
                const vs_type& st = avg.Get_stat();
                stat.insert(stat.end(), st.begin(), st.end());
            }
        }  // end of parallel for (over groups)

        // get average in avF
        for (unsigned tp = 0; tp < T; ++tp) {
            avF[tp] /= ng;
        }
        if (gen < numgen - 1) {
            // if not final generation, get offspring from the current pop, and
            // then put these into pop, forming a new generation by replacing
            // individuals
            vind_type offs = SelectReproduce(pop, mr0);
            // copy individuals in offs to pop, assigning group and individual
            // numbers
            for (unsigned gn = 0; gn < ng; ++gn) { // groups
                for (unsigned i = 0; i < N; ++i) { // inds in group
                    unsigned j = gn*N + i;
                    offs[j].phenotype.gnum = gn;   // set group number
                    offs[j].phenotype.Set_inum(i); // set individual number
                    pop[j] = offs[j]; // replace with offspring
                }
            }
        } // final generation, just leave population, to be saved

        // all set to start next generation
        ++PrBar;
    }
    PrBar.Final();
    timer.Stop();
    timer.Display();
    h5_write_pop(id.h5OutName);
    if (learn_hist) {
        h5_write_hist(id.h5HistName);
    }
}

// return vector of Ntot offspring from the parent population in ppop, with
// individual payoff being proportional to the probability of delivering a
// gamete, and using mutation and recombination parameters from mr
Speci::vind_type Speci::SelectReproduce(const vind_type& ppop, mut_rec_type& mr)
{
    vind_type offspr;
    offspr.reserve(Ntot);
    unsigned np = ppop.size(); // we ought to have np == Ntot
    if (np > 0) {
        // get discrete distribution with parental payoffs as weights, taking
        // into account the costs of beta and u
        v_type wei(np);
        for (unsigned i = 0; i < np; ++i) {
            const phen_type& ph = ppop[i].phenotype;
            wei[i] = ph.payoff;
            if (wei[i] < 0.0) {
                wei[i] = 0.0;
            }
        }
        rand_discr dscr(wei.begin(), wei.end());
        // get offspring
        for (unsigned j = 0; j < Ntot; ++j) {
            // find "mother" for individual to be constructed
            unsigned imat = dscr(mr.eng);
            const ind_type& matind = ppop[imat];
            // copy of "maternal gamete"
            gam_type matgam = matind.genotype.gam;
            // find "father" for individual to be constructed
            unsigned ipat = dscr(mr.eng);
            const ind_type& patind = ppop[ipat];
            // copy of "paternal gamete"
            gam_type patgam = patind.genotype.gam;
            // construct "intermediate" diploid
            dip_type dip(matgam, patgam);
            ind_type indi(dip.GetGamete(mr));
            // append new individual to offspr, constructed with gamete (with
            // mutation and recombination) from dip
            offspr.push_back(indi);
        }
    }
    return offspr;
}

void Speci::h5_read_pop(const std::string& infilename)
{
    // read data and put in pop
    h5R h5(infilename);
    std::vector<v_type> gams(Ntot, v_type(num_loci));
    // read gametes
    h5.read_flt_arr("Gam", gams);
    for (unsigned i = 0; i < Ntot; ++i) {
        gam_type& gam = pop[i].genotype.gam;
        for (unsigned l = 0; l < num_loci; ++l) {
            gam[l] = gams[i][l];
        }
    }
    v_type fval(Ntot);
    // Qini
    h5.read_flt("Qini", fval);
    for (unsigned i = 0; i < Ntot; ++i) {
        pop[i].phenotype.Qini = fval[i];
    }
    // alph
    h5.read_flt("alph", fval);
    for (unsigned i = 0; i < Ntot; ++i) {
        pop[i].phenotype.alph = fval[i];
    }
    // v
    h5.read_flt("v", fval);
    for (unsigned i = 0; i < Ntot; ++i) {
        pop[i].phenotype.v = fval[i];
    }
    // q
    h5.read_flt("q", fval);
    for (unsigned i = 0; i < Ntot; ++i) {
        pop[i].phenotype.q = fval[i];
    }
    // Q0
    h5.read_flt("Q0", fval);
    for (unsigned i = 0; i < Ntot; ++i) {
        pop[i].phenotype.Q0 = fval[i];
    }
    // Q1
    h5.read_flt("Q1", fval);
    for (unsigned i = 0; i < Ntot; ++i) {
        pop[i].phenotype.Q1 = fval[i];
    }
    // p
    h5.read_flt("p", fval);
    for (unsigned i = 0; i < Ntot; ++i) {
        pop[i].phenotype.p = fval[i];
    }
    // payoff
    h5.read_flt("payoff", fval);
    for (unsigned i = 0; i < Ntot; ++i) {
        pop[i].phenotype.payoff = fval[i];
    }
    // R
    h5.read_flt("R", fval);
    for (unsigned i = 0; i < Ntot; ++i) {
        pop[i].phenotype.R = fval[i];
    }
    // std::vector to hold int (or bool) member
    std::vector<int> ival(Ntot);
    // kP
    h5.read_int("kP", ival);
    for (unsigned i = 0; i < Ntot; ++i) {
        pop[i].phenotype.kP = ival[i];
    }
    // std::vector to hold unsigned int member
    ui_type uival(Ntot);
    // u
    h5.read_uint("u", uival);
    for (unsigned i = 0; i < Ntot; ++i) {
        pop[i].phenotype.u = uival[i];
    }
    // nFor
    h5.read_uint("nFor", uival);
    for (unsigned i = 0; i < Ntot; ++i) {
        pop[i].phenotype.nFor = uival[i];
    }
    // nS
    h5.read_uint("nS", uival);
    for (unsigned i = 0; i < Ntot; ++i) {
        pop[i].phenotype.nS = uival[i];
    }
    // inum
    h5.read_uint("inum", uival);
    for (unsigned i = 0; i < Ntot; ++i) {
        pop[i].phenotype.inum = uival[i];
    }
    // gnum
    h5.read_uint("gnum", uival);
    for (unsigned i = 0; i < Ntot; ++i) {
        pop[i].phenotype.gnum = uival[i];
    }
    // female
    h5.read_int("female", ival);
    for (unsigned i = 0; i < Ntot; ++i) {
        pop[i].phenotype.female = ival[i];
    }
    // alive
    h5.read_int("alive", ival);
    for (unsigned i = 0; i < Ntot; ++i) {
        pop[i].phenotype.alive = ival[i];
    }
}

void Speci::h5_write_pop(const std::string& outfilename) const
{
    h5W h5(outfilename);
    std::vector<v_type> gams(Ntot, v_type(num_loci));
    // write gametes
    for (unsigned i = 0; i < Ntot; ++i) {
        const gam_type& gam = pop[i].genotype.gam;
        for (unsigned l = 0; l < num_loci; ++l) {
            gams[i][l] = gam[l];
        }
    }
    h5.write_flt_arr("Gam", gams);
    // write members of phenotypes
    // std::vector to hold flt member
    v_type fval(Ntot);
    // Qini
    std::transform(pop.begin(), pop.end(), fval.begin(),
                   [](const ind_type& i) -> flt
                   { return i.phenotype.Qini; });
    h5.write_flt("Qini", fval);
    // alph
    std::transform(pop.begin(), pop.end(), fval.begin(),
                   [](const ind_type& i) -> flt
                   { return i.phenotype.alph; });
    h5.write_flt("alph", fval);
    // v
    std::transform(pop.begin(), pop.end(), fval.begin(),
                   [](const ind_type& i) -> flt
                   { return i.phenotype.v; });
    h5.write_flt("v", fval);
    // q
    std::transform(pop.begin(), pop.end(), fval.begin(),
                   [](const ind_type& i) -> flt
                   { return i.phenotype.q; });
    h5.write_flt("q", fval);
    // Q0
    std::transform(pop.begin(), pop.end(), fval.begin(),
                   [](const ind_type& i) -> flt
                   { return i.phenotype.Q0; });
    h5.write_flt("Q0", fval);
    // Q1
    std::transform(pop.begin(), pop.end(), fval.begin(),
                   [](const ind_type& i) -> flt
                   { return i.phenotype.Q1; });
    h5.write_flt("Q1", fval);
    // p
    std::transform(pop.begin(), pop.end(), fval.begin(),
                   [](const ind_type& i) -> flt
                   { return i.phenotype.p; });
    h5.write_flt("p", fval);
    // payoff
    std::transform(pop.begin(), pop.end(), fval.begin(),
                   [](const ind_type& i) -> flt
                   { return i.phenotype.payoff; });
    h5.write_flt("payoff", fval);
    // R
    std::transform(pop.begin(), pop.end(), fval.begin(),
                   [](const ind_type& i) -> flt
                   { return i.phenotype.R; });
    h5.write_flt("R", fval);
    // std::vector to hold int (or bool) member
    std::vector<int> ival(Ntot);
    // kP
    std::transform(pop.begin(), pop.end(), ival.begin(),
                   [](const ind_type& i) -> int
                   { return i.phenotype.kP; });
    h5.write_int("kP", ival);
    // std::vector to hold unsigned int member
    ui_type uival(Ntot);
    // u
    std::transform(pop.begin(), pop.end(), uival.begin(),
                   [](const ind_type& i) -> unsigned
                   { return i.phenotype.u; });
    h5.write_uint("u", uival);
    // nFor
    std::transform(pop.begin(), pop.end(), uival.begin(),
                   [](const ind_type& i) -> unsigned
                   { return i.phenotype.nFor; });
    h5.write_uint("nFor", uival);
    // nS
    std::transform(pop.begin(), pop.end(), uival.begin(),
                   [](const ind_type& i) -> unsigned
                   { return i.phenotype.nS; });
    h5.write_uint("nS", uival);
    // inum
    std::transform(pop.begin(), pop.end(), uival.begin(),
                   [](const ind_type& i) -> unsigned
                   { return i.phenotype.inum; });
    h5.write_uint("inum", uival);
    // gnum
    std::transform(pop.begin(), pop.end(), uival.begin(),
                   [](const ind_type& i) -> unsigned
                   { return i.phenotype.gnum; });
    h5.write_uint("gnum", uival);
    // female
    std::transform(pop.begin(), pop.end(), ival.begin(),
                   [](const ind_type& i) -> int
                   { return i.phenotype.female; });
    h5.write_int("female", ival);
    // alive
    std::transform(pop.begin(), pop.end(), ival.begin(),
                   [](const ind_type& i) -> int
                   { return i.phenotype.alive; });
    h5.write_int("alive", ival);
    // write polarisation indices
    h5.write_flt("avF", avF);
}

void Speci::h5_write_hist(const std::string& histfilename) const
{
    h5W h5(histfilename);
    unsigned hlen = stat.size();
    // std::vector to hold unsigned int member
    ui_type uival(hlen);
    // gnum
    std::transform(stat.begin(), stat.end(), uival.begin(),
                   [](const stat_type& st) -> unsigned
                   { return st.gnum; });
    h5.write_uint("gnum", uival);
    // tstep
    std::transform(stat.begin(), stat.end(), uival.begin(),
                   [](const stat_type& st) -> unsigned
                   { return st.tstep; });
    h5.write_uint("tstep", uival);
    // i
    std::transform(stat.begin(), stat.end(), uival.begin(),
                   [](const stat_type& st) -> unsigned
                   { return st.i; });
    h5.write_uint("i", uival);
    // u
    std::transform(stat.begin(), stat.end(), uival.begin(),
                   [](const stat_type& st) -> unsigned
                   { return st.u; });
    h5.write_uint("u", uival);
    // std::vector to hold flt member
    v_type fval(hlen);
    // p
    std::transform(stat.begin(), stat.end(), fval.begin(),
                   [](const stat_type& st) -> flt
                   { return st.p; });
    h5.write_flt("p", fval);
    // R
    std::transform(stat.begin(), stat.end(), fval.begin(),
                   [](const stat_type& st) -> flt
                   { return st.R; });
    h5.write_flt("R", fval);
    // delt
    std::transform(stat.begin(), stat.end(), fval.begin(),
                   [](const stat_type& st) -> flt
                   { return st.delt; });
    h5.write_flt("delt", fval);
    // Q0
    std::transform(stat.begin(), stat.end(), fval.begin(),
                   [](const stat_type& st) -> flt
                   { return st.Q0; });
    h5.write_flt("Q0", fval);
    // Q1
    std::transform(stat.begin(), stat.end(), fval.begin(),
                   [](const stat_type& st) -> flt
                   { return st.Q1; });
    h5.write_flt("Q1", fval);
}
