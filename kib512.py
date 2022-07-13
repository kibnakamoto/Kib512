import numpy as np
import sys

# 64-bit byteswap
def bs64(x):
    return np.uint64((x & np.uint64(0xff00000000000000)) >> np.uint64(56) |
                     (x & np.uint64(0x00ff000000000000)) >> np.uint64(40) |
                     (x & np.uint64(0x0000ff0000000000)) >> np.uint64(24) |
                     (x & np.uint64(0x000000ff00000000)) >> np.uint64(8)  |
                     (x & np.uint64(0x00000000ff000000)) << np.uint64(8)  |
                     (x & np.uint64(0x0000000000ff0000)) << np.uint64(24) |
                     (x & np.uint64(0x000000000000ff00)) << np.uint64(40) |
                     (x & np.uint64(0x00000000000000ff)) << np.uint64(56))

# compression idea: use right rotate and left rotate unieqly so that values that
# are changed do not interfere(the same bit shouldn't be both right rotated
# and left rotated) since that will keep data uncompressed. Instead of
# brute forcing, first generate the values I want to rotate and shift with.
# by using sense and a different program. Find a way to give matrix value even
# if data is equal to zero. Insert a backdoor if possible... 

# unsigned 64-bit bitwise right-rotate
def rr(x: np.uint64(), n: int):
    return (x >> n)|(x << (64-n)) & 0xffffffffffffffff

# unsigned 64-bit bitwise left-rotate
def lr(x: np.uint64(), n: int):
    return (x << n)|(x >> (64-n)) & 0xffffffffffffffff

class Kib512:
    CONST_MATRIX = [
        [0x159a300c24f5dc1e, 0x178f1324ac498bfc, 0x10e55f6926814653,
         0x1963cf80d6276ccc,0x10980aa9695d807a, 0x1703f56bbe1f7f5e,
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
    
    hash_ = {
        0x5287173768046659, 0x1388c1a81885db29, 0xc582055c7b0f1a24,
        0x9be22d058c2ae082, 0xa36b2344c3b2e0d0, 0x8b74c2c08074d4b1,
        0x2a1c87eb1011bd80, 0x64edcba54a4d0a5a
    }
    
    
    """ pre-processing"""
    def prep_kib512(self, inp):
        length = len(inp)
        padlen = (((512-((length*8)+1)-128) % 512)-7)//8
        
        # matrix column height
        self.matrix_colh = (padlen + length + 17)//8;
        
        self.matrix = np.eye(8, self.matrix_colh) # 2-d matrix
        inp+=chr(0x90) # add delimeter after end of data
        inp+='0'*padlen # pad
        inp+=hex(length)[2:].zfill(16) # add length
        arr = np.array(bytearray(inp.encode('charmap')), dtype=np.uint8)
        for i in range(8):
            for j in range(self.matrix_colh):
                self.matrix[i,j] = arr[j+i*self.matrix_colh]
        
        return self.matrix
    
    """ pre-compression """
    def prec_kib512(self):
        self.block_count = self.matrix_colh//8
        self.manip_m = np.zeros((self.block_count, 8, 8), dtype=np.uint64)
        
        mmi = np.uint64(0) # manipulation matrix i
        mmj = np.int32(0) # manipulation matrix j
        for i in range(8):
            for j in range(self.block_count):
                temp = np.uint64(0)
                
                # add matrix values while avoiding repetition of values
                for k in range(8):
                    temp |= np.uint64(self.matrix[i,k+j*8]) << np.uint64(56-k*8)
                
                print(hex(temp)[2:])
                # data in 3-d matrix has to be big-endian
                if(sys.byteorder[0] == 'l'):
                    self.manip_m[mmi,0,mmj] = np.uint64(bs64(temp & 0xffffffffffffffff))
                else:
                    self.manip_m[mmi,0,mmj] = np.uint64(temp & 0xffffffffffffffff)
                mmj = np.int32((mmj+1)%8)
                
            if mmj%8==0:
                mmi+=np.uint64(1)
        
        for i in range(self.block_count):
            for j in range(1,8):
                for k in range(8):
                    # get the first values from the 2d matrix as the value to generate
                    # the row values if you consider that this matrix has stacked
                    # rows and colums. Test until values seem like they have zero
                    # patters and seem completely random
                    
                    # values that is used to generate rest of the matrix
                    temp = self.manip_m[i,j-1,k]
                    self.manip_m[i,j,k]

        for i in range(self.block_count):
            for j in range(1,8):
                for k in range(8):
                    # print(hex(self.manip_m[i,j,k])[2:])
                    pass
        
        return self.manip_m
    
    """ compression function """
    def comp_kib512(self):
        pass
    
    def __call__(self):
        pass

inp = "abcdefghqwertyuioplkjhgfdsazxcvbnm1234567890!@#$%^&*()\\/"
kib512 = Kib512()
kib512.prep_kib512(inp)
kib512.prec_kib512()
