#ifndef HDF5CODE_HPP
#define HDF5CODE_HPP

#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <string>

// The Speci program runs simulations of learning in a social network
// Copyright (C) 2022  Olof Leimar
// See Readme.md for copyright notice

//************************* Class hdf5R **********************************

// This class is used to read a population of individuals from a HDF5 file. The
// file is assumed to have number of datasets in its root group, with names and
// content (data spaces) corresponding to the data fields in the individuals.
// The first dimension of each dataset is assumed to be equal to the number of
// individuals in the population.

class h5R {
public:
    using flt = double;
    using v_type = std::vector<flt>;
    using vv_type = std::vector<v_type>;
    using ui_type = std::vector<unsigned>;
    using i_type = std::vector<int>;
    // using flt_arr = boost::multi_array<flt, 2>;
    // using flt_arr = Array<flt, Dynamic, Dynamic, RowMajor>;
    // constructor opens file for reading
    h5R(std::string in_name) : file(in_name, HighFive::File::ReadOnly) {}
    void read_flt(std::string ds_name, v_type& dat);
    // void read_flt_arr(std::string ds_name, flt_arr& dat);
    // void read_flt_arr(std::string ds_name, flt* dat);
    void read_flt_arr(std::string ds_name, vv_type& dat);
    void read_uint(std::string ds_name, ui_type& dat);
    void read_int(std::string ds_name, i_type& dat);
private:
    // H5::H5File file;  // file (closed when object goes out of scope)
    HighFive::File file;
};

//*************************** Class h5W **********************************

// This class is used to write a population of individuals to a HDF5 file,
// overwriting any content if the file already exists. The data is written as a
// number of datasets in the files root group, with names and content
// (data spaces) corresponding to the data fields in the individuals. The first
// dimension of each dataset is equal to the number of individuals in the
// population.

class h5W {
public:
    using flt = double;
    using v_type = std::vector<flt>;
    using vv_type = std::vector<v_type>;
    using ui_type = std::vector<unsigned>;
    using i_type = std::vector<int>;
    // constructor opens file for writing, truncating any previous file/content
    h5W(std::string out_name) :
        // file(out_name.c_str(), H5F_ACC_TRUNC), num_inds{n_inds} {}
        file(out_name,
             HighFive::File::ReadWrite |
             HighFive::File::Create |
             HighFive::File::Truncate) {}
    void write_flt(std::string ds_name, const v_type& dat);
    void write_flt_arr(std::string ds_name, const vv_type& dat);
    // void write_flt_arr(std::string ds_name, flt** dat,
    //                    std::size_t num_cols);
    void write_uint(std::string ds_name, const ui_type& dat);
    void write_int(std::string ds_name, const i_type& dat);
private:
    // H5::H5File file;        // file (closed when object goes out of scope)
    HighFive::File file;
    // std::size_t num_inds;   // number of individuals in population
};

#endif // HDF5CODE_HPP
