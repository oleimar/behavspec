#ifndef INDIVIDUAL_HPP
#define INDIVIDUAL_HPP

#include <utility>
#include <string>
#include <ostream>
#include <istream>
#include <sstream>

// The Speci program runs simulations of learning in a social network
// Copyright (C) 2022  Olof Leimar
// See Readme.md for copyright notice

//************************ Struct Individual *****************************

// This struct represents an individual with genotype of type GenType,
// phenotype of type PhenType

template<typename GenType, typename PhenType>
struct Individual {
// public:
    using gen_type = GenType;
    using gam_type = typename gen_type::gam_type;
    using mut_rec_type = typename gen_type::mut_rec_type;
    // using rho_vec_type = typename gen_type::rho_vec_type;
    using phen_type = PhenType;
    Individual() {}
    Individual(unsigned g_size, gen_type&& g, phen_type&& ph) :
        genotype{g},
        phenotype{g_size, ph} {}
    // Construct individual from one gamete
    Individual(unsigned g_size, gam_type&& gam) :
        genotype(std::forward<gam_type>(gam)),
        phenotype(g_size, genotype) {}
    Individual(const gam_type& gam) :
        genotype(gam),
        phenotype(genotype) {}
    // Construct individual from maternal and paternal gametes
    Individual(unsigned g_size, gam_type&& mat_gam, gam_type&& pat_gam) :
        genotype(std::forward<gam_type>(mat_gam),
                 std::forward<gam_type>(pat_gam)),
        phenotype(g_size, genotype) {}
    void Assign(gam_type&& gam);
    void Assign(gam_type&& mat_gam, gam_type&& pat_gam);
    void Assign(gen_type&& g, phen_type&& ph);
    gam_type GetGamete(mut_rec_type& mr) const
    { return genotype.GetGamete(mr); }
    // gam_type GetGamete(mut_rec_type& mr, const rho_vec_type& rho) const
    // { return genotype.GetGamete(mr, rho); }
    const gen_type& Genotype() const { return genotype; }
    gen_type& Genotype() { return genotype; }
    const phen_type& Phenotype() const { return phenotype; }
    phen_type& Phenotype() { return phenotype; }
    unsigned SubPopNum() const { return phenotype.spn; }
    void SetSubPopN(unsigned a_spn) { phenotype.spn = a_spn; }
    bool Alive() const { return phenotype.alive; }
    void SetAlive() { phenotype.alive = true; }
    void SetDead() { phenotype.alive = false; }
    bool Female() const { return phenotype.Female(); }
    void SetFemale(bool female) { phenotype.female = female; }
    static std::string ColHeads(unsigned n_loci, unsigned g_size);
    // public data members
    gen_type genotype;
    phen_type phenotype;
};

// Construct individual from one gamete and a spn
template<typename GenType, typename PhenType>
void Individual<GenType, PhenType>::Assign(gam_type&& gam)
{
    genotype.Assign(gam);
    phenotype.Assign(genotype);
}

// Construct individual by assigning maternal and paternal gametes
template<typename GenType, typename PhenType>
void Individual<GenType, PhenType>::Assign(gam_type&& mat_gam,
    gam_type&& pat_gam)
{
    genotype.Assign(mat_gam, pat_gam);
    phenotype.Assign(genotype);
}

// Construct an individual by assigning all its data
template<typename GenType, typename PhenType>
void Individual<GenType, PhenType>::Assign(gen_type&& g, phen_type&& ph)
{
    genotype = g;
    phenotype = ph;
}


template<typename GenType, typename PhenType>
std::string Individual<GenType, PhenType>::ColHeads(unsigned n_loci,
                                                    unsigned g_size)
{
    std::string col_hds = gen_type::ColHeads(n_loci);
    col_hds += "\t";
    col_hds += phen_type::ColHeads(g_size);
    return col_hds;
}


//---------------------------------------------------------------------
// Ouput and input of Individual<T> objects

template<typename GenType, typename PhenType>
std::ostream& operator<<(std::ostream& ostr,
                         const Individual<GenType, PhenType>& ind)
{
    ostr << ind.Genotype() << '\t'
         << ind.Phenotype();
    return ostr;
}

template<typename GenType, typename PhenType>
std::istream& operator>>(std::istream& istr,
                         Individual<GenType, PhenType>& ind)
{
    istr >> ind.Genotype()
         >> ind.Phenotype();
    return istr;
}


#endif // INDIVIDUAL_HPP
