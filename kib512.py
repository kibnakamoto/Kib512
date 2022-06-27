from numpy import *
import qutip as qt

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
        # calculate the first 64 primes to the power of the last 64 primes
        # then right-shift the closest byte size of those values
        const_matrix[i].append((primes[primes_index]**primes[primes_index+ \
                             int(len(primes)/2)]>>int((len(bin(primes[primes_index]** \
                                                          primes[primes_index+ \
                                                          int(len(primes)/2)]) \
                                                         [2:])>>3)))%2**64)
        primes_index+=1; # index of prime array

# for r in const_matrix:
    # print(r)
    # for c in r:
    #     print(hex(c))

# calculate the 8 64-bit values for final hash's starting value
# zx
x_tensordot_z = qt.tensor(qt.sigmax(), qt.sigmaz())
z_tensordot_x = qt.tensor(qt.sigmaz(), qt.sigmax())

# xy
x_tensordot_y = qt.tensor(qt.sigmax(),qt.sigmay())
y_tensordot_x = qt.tensor(qt.sigmay(),qt.sigmax())

# yz
y_tensordot_z = qt.tensor(qt.sigmay(), qt.sigmaz())
z_tensordot_y = qt.tensor(qt.sigmaz(), qt.sigmay())

# zxy
xz_tensordot_y = qt.tensor(x_tensordot_z, qt.sigmay())
zx_tensordot_y = qt.tensor(z_tensordot_x, qt.sigmay())
zy_tensordot_x = qt.tensor(z_tensordot_y, qt.sigmax())
yz_tensordot_x = qt.tensor(y_tensordot_z, qt.sigmax())
xy_tensordot_z = qt.tensor(x_tensordot_y, qt.sigmaz())
yx_tensordot_z = qt.tensor(y_tensordot_x, qt.sigmaz())
xx_tensordot_x = qt.tensor(qt.sigmax(), qt.sigmax(), \
                           qt.sigmax())
zz_tensordot_z = qt.tensor(qt.sigmaz(), qt.sigmaz(), \
                           qt.sigmaz())
yy_tensordot_y = qt.tensor(qt.sigmay(), qt.sigmay(), \
                           qt.sigmay())
var_xxx = ""

for i in xx_tensordot_x:
    for j in i:
        for k in j:
            var_xxx += str(abs(~int(k)%2))
# bitwise and with 129th prime number that went through the same process as primes list. If they aren't all completly unique, use up to the 134th prime number.
print(hex((int(var_xxx,2)&const_matrix[0][0])%2**64))
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

# use tensor(sigmax() sigmaY) from qutip library etc. for more readable output
