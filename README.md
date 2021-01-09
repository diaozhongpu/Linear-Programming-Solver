# Linear-Programming-Solver
Final Project for Applied Operations Research, ZJU 2020-2021 Autumn-Winter

Yanbo Jiang, Zhongpu Diao, and Ziyue Zhan

## Input Specification

```
n m
c_1 ... c_n
a_11 ... a_1n b_1 d_1
......
a_m1 ... a_mn b_m d_m
e_1 ... e_n
```

n : variable number

m : constraint number

c_i  : subjective function

a_ij : coefficient matrix

b_i  : constraint coefficient

d_i  : 

	 = -1, <=

	 =  0, ==

	 =  1, >=

e_i  : x_i constraint

	 = -1, x_i<=0

	 =  0, no constraint for x_i

	 =  1, x_i>=0



## Output Specification

```
k
z
x_1 ... x_n
```


k : solution set type

  = -1, no solution

  =  0, no finite solution

  =  1, with finite solution,

(Optional, when having finite solution)

z   : subjective function optimal

x_i : variable values


## Algorithms Overview



### Simplex Method



### Dual





## Program Module














