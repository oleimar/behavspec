# functions to read data from hdf5 file created by simulation program

# return data table for "one-column" phenotype data
h5_dt <- function(hf_name) {
    require(hdf5r)
    require(data.table)
    f.h5 <- H5File$new(hf_name, mode = "r")
    Qini <- f.h5[["Qini"]][]
    alph <- f.h5[["alph"]][]
    ppf <- f.h5[["ppf"]][]
    q <- f.h5[["q"]][]
    Q0 <- f.h5[["Q0"]][]
    Q1 <- f.h5[["Q1"]][]
    p <- f.h5[["p"]][]
    s <- f.h5[["s"]][]
    pna <- f.h5[["pna"]][]
    payoff <- f.h5[["payoff"]][]
    R <- f.h5[["R"]][]
    kC <- f.h5[["kC"]][]
    u <- f.h5[["u"]][]
    inum <- f.h5[["inum"]][] + 1
    gnum <- f.h5[["gnum"]][] + 1
    female <- f.h5[["female"]][]
    alive <- f.h5[["alive"]][]
    f.h5$close_all()
    data.table(Qini = Qini, alph = alph, ppf = ppf,
               q = q, Q0 = Q0, Q1 = Q1, p = p, payoff = payoff, R = R,
               s = s, pna = pna, kC = kC, u = u,
               inum = inum, gnum = gnum, female = female, alive = alive)
}

# return matrix where each row is a maternal gamete value
h5_mat_gam <- function(hf_name) {
    require(hdf5r)
    f.h5 <- H5File$new(hf_name, mode = "r")
    mat_gam <- t(f.h5[["MatGam"]][,])
    f.h5$close_all()
    mat_gam
}

# return matrix where each row is a paternal gamete value
h5_pat_gam <- function(hf_name) {
    require(hdf5r)
    f.h5 <- H5File$new(hf_name, mode = "r")
    pat_gam <- t(f.h5[["PatGam"]][,])
    f.h5$close_all()
    pat_gam
}

# return vector of average polarisation indices
h5_avF <- function(hf_name) {
    require(hdf5r)
    f.h5 <- H5File$new(hf_name, mode = "r")
    avF <- f.h5[["avF"]][]
    f.h5$close_all()
    avF
}

# return data table for learning history data
h5_hdt <- function(hf_name) {
    require(hdf5r)
    require(data.table)
    f.h5 <- H5File$new(hf_name, mode = "r")
    gnum <- f.h5[["gnum"]][] + 1
    tstep <- f.h5[["tstep"]][] + 1
    i <- f.h5[["i"]][] + 1
    u <- f.h5[["u"]][]
    p <- f.h5[["p"]][]
    R <- f.h5[["R"]][]
    delt <- f.h5[["delt"]][]
    Q0 <- f.h5[["Q0"]][]
    Q1 <- f.h5[["Q1"]][]
    f.h5$close_all()
    data.table(gnum = gnum, tstep = tstep, i = i, u = u,
               p = p, R = R, delt = delt, Q0 = Q0, Q1 = Q1)
}
