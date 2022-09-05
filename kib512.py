import numpy as np
from copy import deepcopy


# compression idea: use right rotate and left rotate unieqly so that values that
# are changed do not interfere(the same bit shouldn't be both right rotated
# and left rotated) since that will keep data uncompressed. Instead of
# brute forcing, first generate the values I want to rotate and shift with.
# by using sense and a different program. Find a way to give matrix value even
# if data is equal to zero. Insert a backdoor if possible... 

# unsigned 64-bit bitwise right-rotate
def rr(x: int, n: int):
    return (x >> n) | (x << (64 - n)) % 0xffffffffffffffc5  # modulo is p of the tckp64k1 curve


# unsigned 64-bit bitwise right-rotate
def lr(x: int, n: int):
    return (x << n)%0xffffffffffffffc5 | (x >> (64 - n))  # modulo is p of the tckp64k1 curve

def point_add(xp, yp, xq, yq, p, a):
    # equation for private key and pointG
    # find lambda
    if yp == yq or xp == xq:
        __lambda = ((3 * (xp ** 2) + a) * pow(2 * yp, -1,
                                              p)) % p
    else:
        __lambda = ((yq - yp) * pow(xq - xp, -1, p)) % p
    xr = (__lambda ** 2 - xp - xq) % p
    yr = (__lambda * (xp - xr) - yp) % p
    return (xr % p, yr % p)

def point_double(x, y, p, a):
    __lambda = ((3 * (x ** 2) + a) * pow(2 * y, -1, p)) % p
    xr = __lambda ** 2 - 2 * x

    # use ys's negative values
    yr = __lambda * (x - xr) - y
    return (xr % p, yr % p)

def montgomery_ladder(pointG, prikey, p, a):
    r0 = list(pointG)
    r1 = point_double(r0[0], r0[1], p, a)
    bits = bin(prikey)[3:]
    for i in bits:
        if i == '0':
            r1 = point_add(r0[0], r0[1], r1[0],
                           r1[1], p, a)
            r0 = point_double(r0[0], r0[1], p, a)
        else:
            r0 = point_add(r0[0], r0[1], r1[0],
                           r1[1], p, a)
            r1 = point_double(r1[0], r1[1], p, a)
    return (r0[0], r0[1])


""" polynomial multiplication """
def poly_mul(a, b, p):
    alength, blength = len(a), len(b)
    lst = [0] * (alength + blength - 1)
    for i in range(alength):
        for j in range(blength):
            lst[i + j] = ((a[i] * b[j])%p + lst[i + j]) % p
    return lst


""" polynomial modulo """
def poly_mod(a, f, p):
    lenf = len(f)
    if lenf < 2:
        raise ValueError("f(x) smaller than 2")

    while len(a) >= lenf:
        if a[-1] != 0:
            for i in range(lenf, 1, -1):
                a[-i] = (a[-i] - a[-1] * f[-i]) % p
        a = a[0:len(a) - 1]
    return a


class Tckp64k1:
    def __init__(self):
        self.a = 0x0000000000000000  # taken from SEC Koblitz curves domain parameters
        self.b = 0x0000000000000007  # taken from SEC Koblitz curves domain parameters

        # p is generated without the use of SEC specifications on how to generate p and q
        # since that is generated randomly, there isn't a need to
        # verified as prime using Fermat's Little Theorem in GF(gf_p) where gf_p
        # denotes the potential largest unsigned 64-bit prime number
        self.p = 0xffffffffffffffc5  # largest unsigned 64-bit prime number

        # calculated with sagemath
        self.n = 0xffffffffffffffc6  # order of curve

        # calculated with sagemath
        # generator point
        self.G = (0x0dc2561d0fc35924, 0xc3790f12017191e9)
        self.h = 0x0000000000000001  # co-factor


# only for square matrix multiplication on GF(p) for same size matrices
class Matrix:
    def __init__(self, m, curve, size):
        self.m = m
        self.curve = curve
        self.res = [[0] * size for i in range(size)]
        self.size = size

    # square matrix multiplication on 2d matrices on GF(p)
    # uses polynomial multiplication modulo x^a + x^-b + 1
    # f(x) inspired from operations from modular square root
    def __mul__(self, m2):
        f = [self.curve.a, -self.curve.b, 1]
        for i in range(self.size):
            for j in range(self.size):
                for k in range(self.size):
                    a = [self.curve.a, int(self.m[k][j]), 1]
                    b = [-self.curve.b, int(m2[i][k]), 1]
                    self.res[i][j] = int(self.res[i][j] + poly_mod(poly_mul(a, b, self.curve.p),
                                                                   f, self.curve.p)[1]) % self.curve.p
        return self


class Kib512:
    """ default class initializer """
    def __init__(self):
        # primes used for rotation and shifting
        # chosen so that one big, one small one rotation is made each time
        self.p = (37, 3, 59, 5)
        self.curve = Tckp64k1()
        self.gf_p = self.curve.p
        self.CONST_M = [
            [0x159a300c24f5dc1e, 0x178f1324ac498bfc, 0x10e55f6926814653,
             0x1963cf80d6276ccc, 0x10980aa9695d807a, 0x1703f56bbe1f7f5e,
             0x16eb052c4abc81fe, 0x1f190bb16583b88f], [0x1d6d15eb483d1a05,
             0x1e2c7d2c722534b7, 0x14dcbefdb2afe52d, 0xf5fbc705a80df745,
             0xe1803d026daf382f, 0x4817a3bcb90fa1de, 0x3b8d9167a61165bd,
             0x400b71315483aecf], [0x7e253caa1049536c, 0x7143a810657de4e0,
             0xba3d16917c91aae3, 0x249b4da86d970a00, 0x290ffdedae0a1c66,
             0x8d490ccfccd45ebb, 0x1852f3c77fae1cb7, 0x8f078d043c3be329],
            [0x568d492cc5140ae4, 0x1cd5bbf83190c484, 0x99d177432a40d119,
             0x1f21998dc8c9145e, 0x560d95c1d691bf3f, 0x57131d9cde1633e5,
             0xcccd035627a53b51, 0x13d653d8c3ae0192], [0x5d1e80c0123206a1,
             0x91b1b3f312ac1be8, 0x19bc3eed09a866a2, 0x29ed711bc0303fb1,
             0x35fa6009ed810bd6, 0x69a90fb1b95d63ae, 0x38bedd43ad759c34,
             0xa787e2edc8133e84], [0x552c119614a1f4e5, 0x68d39af3ec9b831f,
             0x3996ac49f898e965, 0x8d79f80825eb64f1, 0x5babe19fdfe4ec70,
             0x7aa6a18653c9a0cc, 0x510e3b69d44d759a, 0x3486650c8e18c529],
            [0x2b2b4051acc6e9a5, 0x49a9817e3d56ff74, 0xc6fd82e30671aeb5,
             0x20be895e589409fd, 0xd8adec17910a1d4f, 0x18882870301f5a0b,
             0x225cb56188423cf3, 0x4b61cc25258ea5c8], [0x63f68eef219e6c01,
             0x7c29a354966d3f55, 0x1cf4c5f3b98ce4de, 0x1476b731259de5f8,
             0x73642b2ec6c2a44b, 0x407adff58f6e054c, 0x811f7cbc6a7e08c7,
             0x277c0323eb636b90]
        ]

        self.H = [
            0x5287173768046659, 0x1388c1a81885db29, 0xc582055c7b0f1a24,
            0x9be22d058c2ae082, 0xa36b2344c3b2e0d0, 0x8b74c2c08074d4b1,
            0x2a1c87eb1011bd80, 0x64edcba54a4d0a5a
        ]

    """ pre-processing """
    def prep_kib512(self, inp):
        length = len(inp)
        padlen = (((512 - ((length * 8) + 1) - 128) % 512) - 7) // 8

        # matrix column height
        self.matrix_colh = (padlen + length + 17) // 8

        self.matrix = np.eye(8, self.matrix_colh)  # 2-d matrix
        inp += chr(0x90)  # add delimeter after end of data
        inp += '0' * padlen  # pad
        inp += hex(length)[2:].zfill(16)  # add length
        arr = np.array(bytearray(inp.encode('charmap')), dtype=np.uint8)
        for i in range(8):
            for j in range(self.matrix_colh):
                self.matrix[i, j] = arr[j + i * self.matrix_colh]

    """ pre-compression """
    def prec_kib512(self):
        self.block_count = self.matrix_colh // 8
        self.manip_m = np.zeros((self.block_count, 8, 8), dtype=np.uint64)

        mmi = np.uint64(0)  # manipulation matrix i
        mmj = np.int32(0)  # manipulation matrix j
        for i in range(8):
            for j in range(self.block_count):
                temp = np.uint64(0)

                # add matrix values while avoiding repetition of values
                for k in range(8):
                    temp |= np.uint64(self.matrix[i, k + j * 8]) << np.uint64(56 - k * 8)
                self.manip_m[mmi, 0, mmj] = np.uint64(temp & 0xffffffffffffffff)

                mmj = np.int32((mmj + 1) % 8)

            if mmj % 8 == 0:
                mmi += np.uint64(1)

        # modular inverse of p
        inv_p = (0xb3e45306eb3e4507, 0x5555555555555542, 0xcbeea4e1a08ad8c4, 0x666666666666664f)

        for i in range(self.block_count):
            for j in range(1, 8):
                for k in range(8):
                    # get the first values from the 2d matrix as the value to generate
                    # the row values if you consider that this matrix has stacked
                    # rows and columns. Test until values seem like they have zero
                    # patters and seem completely random

                    # matrix indexes are chosen within reason, +7,+6 is for
                    # length of matrix and +1 and +0 is for start of the message.
                    tn4 = deepcopy(int(self.manip_m[(i + self.block_count) % self.block_count][j - 1][(k + 7) % 8]))
                    tn3 = deepcopy(int(self.manip_m[(i + self.block_count) % self.block_count][j - 1][(k + 6) % 8]))
                    tn2 = deepcopy(int(self.manip_m[i][j - 1][(k + 1) % 8]))
                    tn1 = deepcopy(int(self.manip_m[i][j - 1][k]))

                    sigma0 = (rr(tn1, self.p[0]) ^ rr(tn2, self.p[1]) ^ (tn3 << self.p[2]) % self.gf_p |
                              (tn4 << self.p[3]) % self.gf_p) % self.gf_p
                    sigma1 = (rr(tn4, self.p[3]) ^ lr(tn1, self.p[0]) ^ (tn2 << self.p[1]) % self.gf_p |
                              (tn3 << self.p[2]) % self.gf_p) % self.gf_p
                    sigma2 = (rr(tn3, self.p[2]) ^ lr(tn4, self.p[3]) ^ (tn1 << self.p[0]) % self.gf_p |
                              (tn2 << self.p[1]) % self.gf_p) % self.gf_p
                    sigma3 = (rr(tn2, self.p[1]) ^ lr(tn3, self.p[2]) ^ (tn4 << self.p[3]) % self.gf_p |
                              (tn1 << self.p[0]) % self.gf_p) % self.gf_p

                    # use arithmetic addition for non-linearity
                    self.manip_m[i, j, k] = (((sigma0 * inv_p[k % 4]) % self.gf_p) +
                                             sigma1 + sigma2 + sigma3) % self.gf_p

    """ compression function """
    def hash_kib512(self):
        copy = [[0] * 8 for i in range(8)]

        # process on const_m for compression
        # const_m is copied in a non-random re-shuffled way. The shuffling is
        # manip_m[i]'s first and last index modulo 8 being the starting index
        # of const_m, this way, the matrix multiplication is not reversible
        # without knowing manip_m, since reversing a matrix is quite costly
        # in terms of performance, brute-forcing won't be a viable way.
        for i in range(self.block_count):
            new_manip_mi = Matrix(self.manip_m[i], self.curve, 8)
            self.h_copy = deepcopy(self.H) # create copy of hash
            n = int(self.manip_m[i, 0, 0]) % 8  # first value of manip_m[i], starting index
            s = int(self.manip_m[i, 7, 7]) % 8  # last value of manip_m[i], starting index
            for j in range(8):
                n = (n + 1) % 8
                for k in range(8):
                    s = (s + 1) % 8
                    copy[j][k] = deepcopy(self.CONST_M[n][s])
            res = (new_manip_mi * copy).res
            muls = []
            for j in range(8):
                muls.append(montgomery_ladder(self.curve.G, int(res[0][j]), self.curve.p, self.curve.a))

            for j in range(1,8):
                # shift counts
                sc0, sc1, sc2, sc3 = self.p[j % 4], self.p[(j + 1) % 4], self.p[(j + 2) % 4], self.p[(j + 3) % 4]

                for k in range(8):
                    # for every j, update index of p so that every 8 value gets the same
                    # shift count. left-shift hash with 64-p[sc]. Add result to tmp0 for
                    # use addition for non-linearity and xor for linearity
                    # the result[0] isn't used since it has the plaintext, Instead,
                    # ECC multiplication is used. Multiply with generator point of curve tckp64k1
                    # this value is added to copy[j][k] * generator point.
                    tmp0 = (rr(self.H[0], sc0) + ((self.H[1] >> sc0) + int(res[j][k])) % self.gf_p) % self.gf_p
                    tmp1 = (rr(self.H[2], sc1) + (self.H[3] >> sc1)) % self.gf_p
                    tmp1 ^= tmp0
                    tmp2 = (rr(self.H[4], sc2) + (self.H[5] >> sc2)) % self.gf_p
                    tmp2 ^= tmp1
                    tmp3 = (rr(self.H[6], sc3) + (self.H[7] >> sc3)) % self.gf_p
                    tmp3 ^= tmp2
                    # print(hex(tmp0)[2:])
                    # ECC point addition with the first 8 results and multiplication with
                    mul = montgomery_ladder(self.curve.G, self.H[k], self.curve.p, self.curve.a)
                    sigma = point_add(mul[0], mul[1], muls[k][0], muls[k][1], self.curve.p, self.curve.a)

                    # add hash to result and tmp3 which has all tmp[i]. For every 8th of the
                    # hash there are 4 tmp variables xor'ed and added to tmp hash copy
                    # index starts from zero so that [0] isn't replaced with [1] on next line
                    self.H[7] = (self.H[6] + sigma[1]) % self.gf_p
                    self.H[6] = self.H[5]
                    self.H[5] = self.H[4]
                    self.H[4] = (self.H[3] + sigma[0]) % self.gf_p
                    self.H[3] = self.H[2]
                    self.H[2] = self.H[1]
                    self.H[1] = self.H[0]
                    self.H[0] = tmp3
            # add hash copy initialized ub beggining of loop to hash
            for j in range(8):
                self.H[j] = (self.H[j] + self.h_copy[j]) % self.gf_p

    def __call__(self):
        return self.H


# inp = "abcdefghqwertyuioplkjhgfdsazxcvbnm1234567890!@#$%^&*()\\/"
inp = "abc"
kib512 = Kib512()
kib512.prep_kib512(inp)
kib512.prec_kib512()
kib512.hash_kib512()
for i in range(8):
    print(hex(kib512.H[i])[2:].zfill(16), end='')