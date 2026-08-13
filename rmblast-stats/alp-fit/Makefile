ALP = ../ALP_1.98_LIB/cpp

CXX = g++
CXXFLAGS = -O3 -I $(ALP)

alp_fit: alp_fit.cpp $(ALP)/libalp.a
	$(CXX) $(CXXFLAGS) -o $@ alp_fit.cpp -L $(ALP) -lalp

$(ALP)/libalp.a:
	$(MAKE) -C $(ALP)

clean:
	rm -f alp_fit
