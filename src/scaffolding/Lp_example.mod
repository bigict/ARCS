var x1;
var x2;
var e_1_2;
var E_1_2;

minimize z: E_1_2;

s.t. con1 : -x1+ x2 + e_1_2 = 327;
s.t. con2 : e_1_2 <= E_1_2;
s.t. con3 : - e_1_2 <= E_1_2;
