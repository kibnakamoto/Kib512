#include <iostream>
#include <stdint.h>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>

#if !defined(UINT8_MAX)
    using uint8_t = unsigned char
#elif !defined(UINT32_MAX)
    using uint32_t = unsigned int
#elif !defined(UINT64_MAX)
    using uint64_t = unsigned long long
#endif

const uint64_t const_matrix[8][8] = {
    {0x159a300c24f5dc1e, 0x178f1324ac498bfc, 0x10e55f6926814653,
     0x1963cf80d6276ccc,0x10980aa9695d807a, 0x1703f56bbe1f7f5e,
     0x16eb052c4abc81fe, 0x1f190bb16583b88f}, {0x1d6d15eb483d1a05,
     0x1e2c7d2c722534b7, 0x14dcbefdb2afe52d, 0xf5fbc705a80df745,
     0xe1803d026daf382f, 0x4817a3bcb90fa1de, 0x3b8d9167a61165bd,
     0x400b71315483aecf}, {0x7e253caa1049536c, 0x7143a810657de4e0,
     0xba3d16917c91aae3, 0x249b4da86d970a00, 0x290ffdedae0a1c66,
     0x8d490ccfccd45ebb, 0x1852f3c77fae1cb7, 0x8f078d043c3be329},
    {0x568d492cc5140ae4, 0x1cd5bbf83190c484, 0x99d177432a40d119,
     0x1f21998dc8c9145e, 0x560d95c1d691bf3f, 0x57131d9cde1633e5,
     0xcccd035627a53b51, 0x13d653d8c3ae0192}, {0x5d1e80c0123206a1,
     0x91b1b3f312ac1be8, 0x19bc3eed09a866a2, 0x29ed711bc0303fb1,
     0x35fa6009ed810bd6, 0x69a90fb1b95d63ae, 0x38bedd43ad759c34,
     0xa787e2edc8133e84}, {0x552c119614a1f4e5, 0x68d39af3ec9b831f,
     0x3996ac49f898e965, 0x8d79f80825eb64f1, 0x5babe19fdfe4ec70, 
     0x7aa6a18653c9a0cc, 0x510e3b69d44d759a, 0x3486650c8e18c529},
    {0x2b2b4051acc6e9a5, 0x49a9817e3d56ff74, 0xc6fd82e30671aeb5,
     0x20be895e589409fd, 0xd8adec17910a1d4f, 0x18882870301f5a0b,
     0x225cb56188423cf3, 0x4b61cc25258ea5c8}, {0x63f68eef219e6c01,
     0x7c29a354966d3f55, 0x1cf4c5f3b98ce4de, 0x1476b731259de5f8, 
     0x73642b2ec6c2a44b, 0x407adff58f6e054c, 0x811f7cbc6a7e08c7,
     0x277c0323eb636b90}
};

// pre-processing of kib512
uint8_t** kib512_prep(std::string input)
{
    uint8_t** matrix = nullptr;
    std::stringstream ss;
    
    // matrix column height
    uint64_t len = input.length();
    uint32_t pad = (((512-((len*8)+1)-128) % 512)-7)/8;
    uint64_t m_ch = (pad + len + 17)/8;
    matrix = new uint8_t*[8];
    for(int c=0;c<8;c++) matrix[c] = new uint8_t[m_ch];
    
    // add delimeter to denote end of data
    input += 0x90;
    
    // padding
    for(int c=0;c<pad;c++) {
        input+='0';
    }
    
    // add length of input in hex to end of input input
    ss << std::setfill('0') << std::setw(16) << std::hex << (len);
    input+=ss.str();
    
    // 1-d array to 2-d matrix
    for(int r=0;r<8;r++) {
        for(int c=0;c<m_ch;c++) {
            matrix[r][c] = input[c*8+r];
        }
    }
    
    return matrix;
}

int main()
{
    std::string in = "abcda";
    uint8_t** m = kib512_prep(in);
    std::cout << "\n\n" << in.length();
    for(int r=0;r<8;r++) {
        for(int c=0;c<8;c++) {
            // std::cout << m[r][c];
        }
    }
    return 0;
}
