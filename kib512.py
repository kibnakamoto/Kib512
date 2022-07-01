import calc_consts_kib512 as consts
import numpy as np

const_matrix = [
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

""" pre-processing"""
# user input as np.uint8 matrix

def prep_kib512(inp):

    length = len(inp)
    padlen = (((512-((len*8)+1)-128) % 512)-7)/8
    
    # matrix column height
    matrix_colh = (padlen + length + 17)/8;
    matrix = np.eye(8,matrix_colh)
    
    inp+='0'*pad # pad
    inp+=hex(inp)[2:].zfill(16) # add length
    
    arr = np.array(bytearray(inp.encode('utf-8')), dtype=np.uint8)
    for i in range(0, 8):
        for j in range(0, matrix_colh):
            matrix[i,j] = arr[i+j*8]
    

""" pre-compression """
# top and bottom row and second top and second below of the input matrix

# endiennes of message

""" compression function """
# how many rounds for compression?

# matrix multipication with bitmasking

# have mutiple mini compression functions for the main one

inp = "abcd"
prep_kib512(inp)
