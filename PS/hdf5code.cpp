#include "hdf5code.hpp"
#include <iostream>
#include <string>

// The Speci program runs simulations of learning in a social network
// Copyright (C) 2022  Olof Leimar
// See Readme.md for copyright notice

//*************************** Class h5R **********************************

void h5R::read_flt(std::string ds_name, v_type& dat)
{
    HighFive::DataSet ds = file.getDataSet(ds_name);
    ds.read(dat);
    // H5::DataSet ds = file.openDataSet(ds_name.c_str());
    // ds.read(dat, H5::PredType::NATIVE_FLOAT);
    // ds.read(dat, H5::PredType::NATIVE_DOUBLE);
}

void h5R::read_flt_arr(std::string ds_name, vv_type& dat)
{
    HighFive::DataSet ds = file.getDataSet(ds_name);
    ds.read(dat);
}

// void h5R::read_flt_arr(std::string ds_name, flt* dat)
// {
//     HighFive::DataSet ds = file.getDataSet(ds_name);
//     ds.read(dat);
//     // H5::DataSet ds = file.openDataSet(ds_name.c_str());
//     // ds.read(dat, H5::PredType::NATIVE_FLOAT);
//     // ds.read(dat, H5::PredType::NATIVE_DOUBLE);
// }

void h5R::read_uint(std::string ds_name, ui_type& dat)
{
    HighFive::DataSet ds = file.getDataSet(ds_name);
    ds.read(dat);
//     H5::DataSet ds = file.openDataSet(ds_name.c_str());
//     ds.read(dat, H5::PredType::NATIVE_UINT);
}

void h5R::read_int(std::string ds_name, i_type& dat)
{
    HighFive::DataSet ds = file.getDataSet(ds_name);
    ds.read(dat);
//     H5::DataSet ds = file.openDataSet(ds_name.c_str());
//     ds.read(dat, H5::PredType::NATIVE_INT);
}


//*************************** Class h5W **********************************

void h5W::write_flt(std::string ds_name, const v_type& dat)
{
    HighFive::DataSet ds =
        file.createDataSet<flt>(ds_name, HighFive::DataSpace::From(dat));
    ds.write(dat);
    // hsize_t dims[1];             // dataset dimensions
    // dims[0] = num_inds;
    // H5::DataSpace dsp(1, dims);  // rank is 1
    // H5::IntType dtype(H5::PredType::NATIVE_FLOAT);
    // // H5::IntType dtype(H5::PredType::NATIVE_DOUBLE);
    // // dtype.setOrder(H5T_ORDER_LE);        // perhaps not needed
    // H5::DataSet ds(file.createDataSet(ds_name.c_str(), dtype, dsp));
    // ds.write(dat, H5::PredType::NATIVE_FLOAT);
    // // ds.write(dat, H5::PredType::NATIVE_DOUBLE);
}

void h5W::write_flt_arr(std::string ds_name, const vv_type& dat)
{
    HighFive::DataSet ds =
        file.createDataSet<flt>(ds_name, HighFive::DataSpace::From(dat));
    ds.write(dat);
}

// void h5W::write_flt_arr(std::string ds_name,
//                         flt** dat,
//                         std::size_t num_cols)
// {
//     std::vector<std::size_t> dims(2);             // dataset dimensions
//     dims[0] = num_inds;
//     dims[1] = num_cols;
//     HighFive::DataSet ds =
//         file.createDataSet<flt>(ds_name, HighFive::DataSpace(dims));
//     ds.write(dat);
//     hsize_t dims[2];             // dataset dimensions
//     H5::DataSpace dsp(2, dims);  // rank is 2
//     H5::IntType dtype(H5::PredType::NATIVE_FLOAT);
//     // H5::IntType dtype(H5::PredType::NATIVE_DOUBLE);
//     // dtype.setOrder(H5T_ORDER_LE);        // perhaps not needed
//     H5::DataSet ds(file.createDataSet(ds_name.c_str(), dtype, dsp));
//     ds.write(dat, H5::PredType::NATIVE_FLOAT);
//     // ds.write(dat, H5::PredType::NATIVE_DOUBLE);
// }

void h5W::write_uint(std::string ds_name, const ui_type& dat)
{
    HighFive::DataSet ds =
        file.createDataSet<unsigned>(ds_name, HighFive::DataSpace::From(dat));
    ds.write(dat);
//     hsize_t dims[1];             // dataset dimensions
//     dims[0] = num_inds;
//     H5::DataSpace dsp(1, dims);  // rank is 1
//     H5::IntType dtype(H5::PredType::NATIVE_UINT);
//     // dtype.setOrder(H5T_ORDER_LE);        // perhaps not needed
//     H5::DataSet ds(file.createDataSet(ds_name.c_str(), dtype, dsp));
//     ds.write(dat, H5::PredType::NATIVE_UINT);
}

void h5W::write_int(std::string ds_name, const i_type& dat)
{
    HighFive::DataSet ds =
        file.createDataSet<int>(ds_name, HighFive::DataSpace::From(dat));
    ds.write(dat);
//     hsize_t dims[1];             // dataset dimensions
//     dims[0] = num_inds;
//     H5::DataSpace dsp(1, dims);  // rank is 1
//     H5::IntType dtype(H5::PredType::NATIVE_INT);
//     // dtype.setOrder(H5T_ORDER_LE);        // perhaps not needed
//     H5::DataSet ds(file.createDataSet(ds_name.c_str(), dtype, dsp));
//     ds.write(dat, H5::PredType::NATIVE_INT);
}
