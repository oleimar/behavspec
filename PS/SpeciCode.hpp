#ifndef SPECICODE_HPP
#define SPECICODE_HPP

// NOTE: for easier debugging one can comment out the following
#ifdef _OPENMP
#define PARA_RUN
#endif

#include "Genotype.hpp"
#include "Phenotype.hpp"
#include "Individual.hpp"
#include "AVgroup.hpp"
#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <random>
#include <algorithm>

// The Speci program runs simulations of learning in a social network
// Copyright (C) 2022  Olof Leimar
// See Readme.md for copyright notice

// An individual has 3 genetically determined traits, Qini, alph, v (see
// Phenotype.hpp), and there is one locus for each trait

//************************* Class SpeciInpData ***************************

// This class is used to 'package' input data in a single place; the
// constructor extracts data from an input file

class SpeciInpData {
public:
    using flt = double;         // the TOML class requires doubles
    using v_type = std::vector<flt>;
    std::size_t max_num_thrds;  // Max number of threads to use
    unsigned num_loci;          // Number of loci in individual's genotype
    unsigned ng;                // Number of groups
    unsigned N;                 // Group size
    unsigned KR;                // Neighbours to the right; degree is 2*KR
    unsigned nScr;              // Maximum number of scroungers at food
    unsigned T;                 // Number of rounds in group interaction
    unsigned T0hist;            // Start round for history
    unsigned numgen;            // Number of generation to simulate
    flt prwr;                   // Probability to "rewire edge" in graph
    flt V0;                     // Baseline payoff per generation
    flt V1;                     // Payoff parameter
    flt V2n;                    // Payoff parameter
    flt V2x;                    // Payoff parameter
    flt qn;                     // Minimum q
    flt qx;                     // Maximum q
    flt prs;                    // Producer success parameter
    flt pmax;                   // Maximum value for probability to use action
    flt bet;                    // parameter for soft-max function
    v_type mut_rate;            // Probability of mutation at each locus
    v_type SD;                  // SD of mutational increments at each locus
    v_type max_val;             // Maximal allelic value at each locus
    v_type min_val;             // Minimal allelic value at each locus
    v_type rho;                 // Recombination rates
    v_type all0;                // Starting allelic values (if not from file)
    bool learn_hist;            // Whether to compute and save learning history
    bool read_from_file;        // Whether to read population from file
    std::string h5InName;       // File name for input of population
    std::string h5OutName;      // File name for output of population
    std::string h5HistName;     // File name for output of learning history

    std::string InpName;  // Name of input data file
    bool OK;              // Whether input data has been successfully read

    SpeciInpData(const char* filename);
};


//***************************** Class Speci ******************************

class Speci {
public:
    // types needed to define individual
    // using mut_rec_type = MutRec<MutIncrBiExp<>>;
    using mut_rec_type = MutRec<MutIncrNorm<>>;
    using gam_type = Gamete<mut_rec_type>;
    using gen_type = Haplotype<gam_type>;
    using dip_type = Diplotype<gam_type>;
    using phen_type = Phenotype<gen_type>;
    using ind_type = Individual<gen_type, phen_type>;
    using stat_type = LearnStat<phen_type>;
    // use std::vector containers for (sub)populations
    using vind_type = std::vector<ind_type>;
    using vph_type = std::vector<phen_type>;
    using avg_type = ActValGroup<phen_type>;
    using flt = double;
    using v_type = std::vector<flt>;
    using ui_type = std::vector<unsigned>;
    using vs_type = std::vector<stat_type>;
    using rand_eng = std::mt19937;
    using vre_type = std::vector<rand_eng>;
    // using rand_int = std::uniform_int_distribution<int>;
    using rand_uni = std::uniform_real_distribution<flt>;
    using rand_norm = std::normal_distribution<flt>;
    using rand_discr = std::discrete_distribution<int>;
    Speci(const SpeciInpData& sid);
    void Run();
    void h5_read_pop(const std::string& infilename);
    void h5_write_pop(const std::string& outfilename) const;
    void h5_write_hist(const std::string& histfilename) const;
private:
    vind_type SelectReproduce(const vind_type& ppop, mut_rec_type& mr);

    SpeciInpData id;
    unsigned num_loci;
    unsigned ng;
    unsigned N;
    unsigned KR;
    unsigned nScr;
    unsigned Ntot;
    unsigned T;
    unsigned T0hist;
    unsigned numgen;
    flt prwr;
    flt V0;
    flt V1;
    flt V2;
    flt V2n;
    flt V2x;
    flt qn;
    flt qx;
    flt prs;
    flt pmax;
    flt bet;
    bool learn_hist;
    std::size_t num_thrds;
    ui_type sds;
    vre_type vre;
    vind_type pop;
    v_type avF;
    vs_type stat;
};

#endif // SPECICODE_HPP
