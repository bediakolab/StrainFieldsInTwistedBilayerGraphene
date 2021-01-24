% vector_equivalence_testing.m

syms a

symflag = 1;
[ dict_function, vsym ] = getLatticePlaneNormalVectors(symflag);
v1sym = [ 0, (3^(1/2)*a)/3];
v2sym = [ a/2, (3^(1/2)*a)/6];

w1sym = v1sym+v2sym;
w2sym = -v1sym+2*v2sym;

k = [1 1];
delta0 = [0.9,-0.3];

deltaf = k(1)*w1sym + k(2)*w2sym - delta0;
