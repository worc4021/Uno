#  NLP written by GAMS Convert at 07/09/02 21:32:07
#  
#  Equation counts
#     Total       E       G       L       N       X       C
#         2       2       0       0       0       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#         4       4       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#         7       4       3       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x2 := 0.51802, >= 1E-6, <= 1;
var x3 := 0.0511, >= 1E-6, <= 1;
var x4 := 0.43088, >= 1E-6, <= 1;

minimize obj: (26.9071667605344*x2 + 41.7710875549227*x3 + 6.30931398488382*x4)
              *log(3.9235*x2 + 6.0909*x3 + 0.92*x4) + 0.668686155614739*x2 - 
              1.14374230885457*x3 + 2.8906196099828*x4 + 9.58716676053442*x2*
              log(x2) + 16.9310875549227*x3*log(x3) + 0.309313984883821*x4*log(
              x4) - 9.58716676053442*x2*log(3.9235*x2 + 6.0909*x3 + 0.92*x4) - 
              16.9310875549227*x3*log(3.9235*x2 + 6.0909*x3 + 0.92*x4) - 
              0.309313984883821*x4*log(3.9235*x2 + 6.0909*x3 + 0.92*x4) + 18.32
              *x2*log(x2) + 25.84*x3*log(x3) + 7*x4*log(x4) - 18.32*x2*log(
              3.664*x2 + 5.168*x3 + 1.4*x4) - 25.84*x3*log(3.664*x2 + 5.168*x3
               + 1.4*x4) - 7*x4*log(3.664*x2 + 5.168*x3 + 1.4*x4) + (4.0643*x2
               + 5.7409*x3 + 1.6741*x4)*log(4.0643*x2 + 5.7409*x3 + 1.6741*x4)
               + 4.0643*x2*log(x2) + 5.7409*x3*log(x3) + 1.6741*x4*log(x4) - 
              4.0643*x2*log(4.0643*x2 + 3.22644664511275*x3 + 1.44980651607875*
              x4) - 5.7409*x3*log(5.31147575751424*x2 + 5.7409*x3 + 
              0.00729924451284409*x4) - 1.6741*x4*log(2.25846661774355*x2 + 
              3.70876916588753*x3 + 1.6741*x4) - 30.9714667605344*x2*log(x2) - 
              47.5119875549227*x3*log(x3) - 7.98341398488382*x4*log(x4);

subject to

e2:    x2 + x3 + x4 = 1;