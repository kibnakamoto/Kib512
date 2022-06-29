from numpy import *
import sys

# find the first 128 prime numbers
primes = []
const_matrix = [[] for i in range(8)] # 4096-bits

inc = 2
while len(primes) != 128:
    if (2**(inc-1)-1) % inc == 0:
        primes.append(inc)
    inc+=1

primes_index = 0;

# calculate the constant values for the iteration matrix
for i in range(8):
    for j in range(8):
        # calculate the first 64 primes to the power of the next 64 primes
        # then right-shift the closest byte size of those values
        temp = (primes[primes_index]** \
                                primes[primes_index+ \
                                int(len(primes)/2)]>> \
                                int((len(bin(primes[primes_index]** \
                                             primes[primes_index+ \
                                             int(len(primes)/2)]) \
                                         [2:])>>3)))
        
        temp = int(hex(temp)[2:18],16)
        
        const_matrix[i].append(temp)
        primes_index+=1; # index of prime array
# 277, 811, 407, 736
# for r in const_matrix:
    # print(r)
    # for c in r:
    #     print(hex(c))

# calculate the 8 64-bit values for final hash's starting value
# zx
sigmax = [[0,1],[1,0]]
sigmay = [[0,-1],[1,0]]
sigmaz = [[1,0],[0,-1]]
x_tensordot_z = kron(sigmax, sigmaz)
z_tensordot_x = kron(sigmaz, sigmax)

# xy
x_tensordot_y = kron(sigmax,sigmay)
y_tensordot_x = kron(sigmay,sigmax)

# yz
y_tensordot_z = kron(sigmay, sigmaz)
z_tensordot_y = kron(sigmaz, sigmay)

# zxy
xz_tensordot_y = kron(x_tensordot_z, sigmay)
zx_tensordot_y = kron(z_tensordot_x, sigmay)
zy_tensordot_x = kron(z_tensordot_y, sigmax)
yz_tensordot_x = kron(y_tensordot_z, sigmax)
xy_tensordot_z = kron(x_tensordot_y, sigmaz)
yx_tensordot_z = kron(y_tensordot_x, sigmaz)
xx_tensordot_x = kron(kron(sigmax, sigmax), \
                           sigmax)
zz_tensordot_z = kron(kron(sigmaz, sigmaz), \
                           sigmaz)
yy_tensordot_y = kron(kron(sigmay, sigmay), \
                           sigmay) # do not include 
var_xxx = ""

# calculate 129th to 137th prime numbers
temp_primes = []
inc = primes[len(primes)-1]+1
while len(temp_primes) != 16: # generate 16 prime numbers
    if (2**(inc-1)-1) % inc == 0:
        temp_primes.append(inc)
    inc+=1

def matrix_to_int(matrix):
    bits = 0
    for i in matrix:
        for j in i:
            bits<<=1
            bits|=abs(~int(j)%2)
    return bits

H = [None]*8 # declare Hash list

# transform temp primes
trnsfm_tmp_ps = []
for i in range(8):
   trnsfm_tmp_ps.append((temp_primes[i]** \
                       temp_primes[i+int(len(temp_primes)/2)]>> \
                       int((len(bin(temp_primes[i]**temp_primes[i+ \
                                    int(len(temp_primes)/2)]) \
                                [2:])>>3)))%2**64)
    
# try matmul matrix_to_int parameter with trnsfm_tmp_ps
H[0] = (matrix_to_int(matmul(xz_tensordot_y,yy_tensordot_y))&
        trnsfm_tmp_ps[0])%2**64
H[1] = (matrix_to_int(matmul(zx_tensordot_y,yy_tensordot_y))&
        trnsfm_tmp_ps[1])%2**64
H[2] = (matrix_to_int(matmul(zy_tensordot_x,yy_tensordot_y))&
        trnsfm_tmp_ps[2])%2**64
H[3] = (matrix_to_int(matmul(yz_tensordot_x,yy_tensordot_y))&
        trnsfm_tmp_ps[3])%2**64
H[4] = (matrix_to_int(matmul(xy_tensordot_z,yy_tensordot_y))&
        trnsfm_tmp_ps[4])%2**64
H[5] = (matrix_to_int(matmul(yx_tensordot_z,yy_tensordot_y))&
        trnsfm_tmp_ps[5])%2**64
H[6] = (matrix_to_int(matmul(xx_tensordot_x,
                             yy_tensordot_y))&
        trnsfm_tmp_ps[6])%2**64
H[7] = (matrix_to_int(matmul(zz_tensordot_z,
                             yy_tensordot_y))&
        trnsfm_tmp_ps[7])%2**64

for i in const_matrix:
    for j in i:
        print(hex(j))

# tensor product of pauli x and pauli z
# 0  0     1  0
# 0  0     0 -1
# 1  0     0  0
# 0 -1     0  0
# tensor product of pauli z and pauli x
# 0  1     0  0
# 1  0     0  0
# 0  0     0 -1
# 0  0    -1  0

# tensor product of pauli x and pauli y
# 0  0     0 -1
# 0  0     1  0
# 0 -1     0  0
# 0  1     0  0
# tensor product of pauli y and pauli x
# 0  0     0 -1
# 0  0    -1  0
# 0  1     0  0
# 1  0     0  0

# tensor product of pauli y and pauli z
# 0  0    -1  0
# 0  0     0  1
# 1  0     0  0
# 0 -1     0  0
# tensor product of pauli z and pauli y
# 0 -1     0  0
# 1  0     0  0
# 0  0     0  1
# 0  0    -1  0

# xz_tensordot_y
# 0  0     0  0     0 -1     0  0
# 0  0     0  0     1  0     0  0
# 0  0     0  0     0  0     0  1
# 0  0     0  0     0  0    -1  0
# 0 -1     0  0     0  0     0  0
# 1  0     0  0     0  0     0  0
# 0  0     0  1     0  0     0  0
# 0  0    -1  0     0  0     0  0

# zx_tensordot_y
# 0  0     0 -1     0  0     0  0
# 0  0     1  0     0  0     0  0
# 0 -1     0  0     0  0     0  0
# 1  0     0  0     0  0     0  0
# 0  0     0  0     0  0     0  1
# 0  0     0  0     0  0    -1  0
# 0  0     0  0     0  1     0  0
# 0  0     0  0    -1  0     0  0

# zy_tensordot_x
# 0  0     0 -1     0  0     0  0
# 0  0    -1  0     0  0     0  0
# 0  1     0  0     0  0     0  0
# 1  0     0  0     0  0     0  0
# 0  0     0  0     0  0     0  1
# 0  0     0  0     0  0     1  0
# 0  0     0  0     0 -1     0  0
# 0  0     0  0    -1  0     0  0

# yz_tensordot_x
# 0  0     0  0     0 -1     0  0
# 0  0     0  0    -1  0     0  0
# 0  0     0  0     0  0     0  1
# 0  0     0  0     0  0     1  0
# 0  1     0  0     0  0     0  0
# 1  0     0  0     0  0     0  0
# 0  0     0 -1     0  0     0  0
# 0  0    -1  0     0  0     0  0

# xy_tensordot_z
# 0  0     0  0     0  0    -1  0
# 0  0     0  0     0  0     0  1
# 0  0     0  0     1  0     0  0
# 0  0     0  0     0 -1     0  0
# 0  0    -1  0     0  0     0  0
# 0  0     0  1     0  0     0  0
# 1  0     0  0     0  0     0  0
# 0 -1     0  0     0  0     0  0

# yx_tensordot_z
# 0  0     0  0     0  0    -1  0
# 0  0     0  0     0  0     0  1
# 0  0     0  0    -1  0     0  0
# 0  0     0  0     0  1     0  0
# 0  0     1  0     0  0     0  0
# 0  0     0 -1     0  0     0  0
# 1  0     0  0     0  0     0  0
# 0 -1     0  0     0  0     0  0

# (xx_tensordot_x bitor yy_tensordot_y bitor xx_tensordot_x) % 0x2

# (zz_tensordot_z bitor yy_tensordot_y bitor xx_tensordot_x) % 0x2

# all the matrices have to be a combination of each other

# 8 final matrices in total, use bitwise not since most numbers are zeros.
# use abs function to convert all to positive integers. perceive the matrix values
# as bits and convert to 8 64-bit unsigned integers and use a similiar logic as 
# how the sha512 const H vector.
