#!/usr/bin/python
from apop import *
from pdb import set_trace

def test_percentiles():
    v = gsl_vector(307)
    for i in range(307):
          v.set(i, i)
    pcts_up    = v.percentiles('u')
    pcts_down  = v.percentiles('d')
    pcts_avg   = v.percentiles('a')
    for i in range(101):
        assert pcts_up[i] >= pcts_down[i]
        assert pcts_up[i] >= pcts_avg[i]
        assert pcts_avg[i] >= pcts_down[i]
    assert (pcts_up[100] == pcts_down[100]) & (pcts_avg[100] == pcts_down[100])
    assert (pcts_up[0] == pcts_down[0]) & (pcts_avg[0] == pcts_down[0])
    assert (pcts_avg[50] == (pcts_down[50] + pcts_up[50])/2)

def test_listwise_delete():
    t1 = apop_data_calloc(0, 10,10)
    t1c = t1.listwise_delete()
    assert t1c.matrix.size1==10
    assert t1c.matrix.size2==10
    t1.vector = gsl_vector_calloc(10)
    t1.set(4,-1, NaN)
    t2c = t1.listwise_delete()
    assert t2c.matrix.size1==9
    t1.set(4,-1, NaN)
    t1.set(7,-1, NaN)
    t3c = t1.listwise_delete()
    assert t3c.matrix.size1==7
    vive = apop_col(t1, 7)
    vive.set_all(NaN)
    assert t1.listwise_delete() == None

def test_skew_and_kurt():
    r  = apop_rng_alloc(38)
    apop_table_exists("t", 1)
    apop_query("create table t(vals)")
    for i in range(1000):
        apop_query("insert into t values(%g)", gsl_rng_uniform(r))
    v    = apop_query_to_vector("select * from t")
    assert (abs(v.var() -apop_query_to_float("select var(vals) from t"))<1e-6)
    assert (abs(v.skew() -apop_query_to_float("select skew(vals) from t"))<1e-6)
    assert (abs(v.kurt() -apop_query_to_float("select kurt(vals) from t"))<1e-5)
    apop_table_exists("t",1)

def test_nan_data():
    apop_text_to_db("test_data_nans", "nandata", 0,1, None)
    avar.apop_opts.db_name_column = "head"
    avar.apop_opts.db_nan =  "\\."
    d  = apop_query_to_data("select * from nandata")
##avar.apop_opts.output_type ='d'
#    apop_data_print(d, "nantest", 0, 0, 's')
#    apop_text_to_db("nantest", "nantest", 0,1, None)
#    d2  = apop_query_to_data("select * from nantest")
    d2  = apop_query_to_data("select * from nandata")
    assert(isNaN(d2.get("second", "c")))
    assert(isNaN(d2.get("third", "b")))
    assert(not d2.get("fourth", "b"))
    del(d2)
    avar.apop_opts.db_nan == "XX"

def wmt(v, v2, w, av, av2, mean):
    assert(v.mean() == v.weighted_mean(None))
    assert(av.mean() == v.weighted_mean(w))
    assert(v.weighted_mean(w) == mean)
    assert(abs(v.var() - v.weighted_var(None))<1e-5)
    assert(abs(v.cov(v2) - v.weighted_cov(v2,None))<1e-5)
    assert(abs(av.var() - v.weighted_var(w))<1e-5)
    assert(abs(av2.cov() - v.weighted_cov(v2,w))<1e-5)
    assert(abs(av.skew_pop() - v.weighted_skew(w))<1e-5)
    assert(abs(av.kurtosis_pop() - v.weighted_kurt(w))<1e-5)

def test_weigted_moments():
  data      = [1,2,3]
  alldata   = [1,2,3]
  data3     = [3,2,1]
  alldata3  = [3,2,1]
  weights   = [1,1,1]
  v          = apop_array_to_vector(data, 3)
  v2         = apop_array_to_vector(data3, 3)
  w          = apop_array_to_vector(weights, 3)
  av         = apop_array_to_vector(alldata, 3)
  av2        = apop_array_to_vector(alldata3, 3)
  wmt(v,v2,w,av,av2,2)
  data2       = [0,1,2,3,4]
  alldata2    = [0,0,0,0,1,1,1,2,2,3]
  data4       = [0,1,3,2,4]
  alldata4    = [0,0,0,0,1,1,1,3,3,2]
  weights2    = [4,3,2,1,0]
  v             = apop_array_to_vector(data2, 5)
  v2            = apop_array_to_vector(data4, 5)
  w             = apop_array_to_vector(weights2, 5)
  av            = apop_array_to_vector(alldata2, 10)
  av2           = apop_array_to_vector(alldata4, 10)
  wmt(v,v2,w,av,av2,1)

TOL= 1e-15
TOL2= 1e-5
TOL3= 1e-9
TOL4= 1e-4

def pontius():
    apop_text_to_db("pontius.dat","d", 0,1, None)
    d        = apop_query_to_data("select y, x, pow(x,2) as p from d")
    est  =  apop_estimate(d, avar.apop_ols)

    assert(abs(est.parameters.get(0, -1) - 0.673565789473684E-03) < TOL)
    assert(abs(est.parameters.get( 1, -1) - 0.732059160401003E-06) < TOL)
    assert(abs(est.parameters.get( 2, -1) - -0.316081871345029E-14)    < TOL)
    cov = est.parameters.more
    assert(abs(cov.get(0, 0) - pow(0.107938612033077E-03,2))    < TOL2)
    assert(abs(cov.get(1, 1) - pow(0.157817399981659E-09,2))    < TOL2)
    assert(abs(cov.get(2, 2) - pow(0.486652849992036E-16,2))    < TOL2)
    cc   = apop_estimate_coefficient_of_determination(est)
    assert(abs(cc.get( "R.sq.*", 0) - 0.999999900178537)    < TOL)
    assert(abs(cc.get( "SSR", 0) - 15.6040343244198)    < TOL3)

def wampler1():
    apop_text_to_db("wampler1.dat","w1", 0,1, None)
    d    = apop_query_to_data("select y, x, pow(x,2) as p2, \
                                pow(x,3) as p3, pow(x,4) as p4, pow(x,5) as p5 from w1")
    est  =  apop_estimate(d, avar.apop_ols)
    for i in range(6):
        assert(abs(est.parameters.get(i, -1) - 1) < TOL4)
    cov = est.parameters.more
    for i in range(6):
        assert(abs(cov.get(i, i)) < TOL2)
    cc   = apop_estimate_coefficient_of_determination(est)
    assert(abs(cc.get( "R.sq.*", 0) - 1)    < TOL)

def numacc4():
    d  = apop_text_to_data("numacc4.dat", 0, 0)
    v = apop_col(d, 0)
    assert(v.mean() == 10000000.2)
    assert(v.var()*(v.size -1)/v.size - 0.01 < TOL3)

print "Testing Python interface",
test_skew_and_kurt()
test_percentiles()
pontius()
wampler1()
numacc4()
test_listwise_delete()
test_nan_data()
print "Done."
