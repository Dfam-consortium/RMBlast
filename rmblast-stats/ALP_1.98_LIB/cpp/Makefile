CXXFLAGS = -O3

OBJ = njn_dynprogprob.o njn_dynprogproblim.o njn_dynprogprobproto.o	\
njn_ioutil.o njn_localmaxstat.o njn_localmaxstatmatrix.o		\
njn_localmaxstatutil.o njn_random.o sls_alignment_evaluer.o sls_alp.o	\
sls_alp_data.o sls_alp_regression.o sls_alp_sim.o sls_basic.o		\
sls_pvalues.o

libalp.a: ${OBJ}
	ar -r $@ $?

clean:
	rm -f *.o *.a
