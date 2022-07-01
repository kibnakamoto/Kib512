#include <iostream>
#include <stdint.h>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <bit>

#if !defined(UINT8_MAX)
    using uint8_t = unsigned char
#elif !defined(UINT32_MAX)
    using uint32_t = unsigned int
#elif !defined(UINT64_MAX)
    using uint64_t = unsigned long long
#endif

const uint64_t const_matrix[8][8] = {
    {0x159a300c24f5dc1eULL, 0x178f1324ac498bfcULL, 0x10e55f6926814653ULL,
     0x1963cf80d6276cccULL,0x10980aa9695d807aULL, 0x1703f56bbe1f7f5eULL,
     0x16eb052c4abc81feULL, 0x1f190bb16583b88fULL}, {0x1d6d15eb483d1a05ULL,
     0x1e2c7d2c722534b7ULL, 0x14dcbefdb2afe52dULL, 0xf5fbc705a80df745ULL,
     0xe1803d026daf382fULL, 0x4817a3bcb90fa1deULL, 0x3b8d9167a61165bdULL,
     0x400b71315483aecfULL}, {0x7e253caa1049536cULL, 0x7143a810657de4e0ULL,
     0xba3d16917c91aae3ULL, 0x249b4da86d970a00ULL, 0x290ffdedae0a1c66ULL,
     0x8d490ccfccd45ebbULL, 0x1852f3c77fae1cb7ULL, 0x8f078d043c3be329ULL},
    {0x568d492cc5140ae4ULL, 0x1cd5bbf83190c484ULL, 0x99d177432a40d119ULL,
     0x1f21998dc8c9145eULL, 0x560d95c1d691bf3fULL, 0x57131d9cde1633e5ULL,
     0xcccd035627a53b51ULL, 0x13d653d8c3ae0192ULL}, {0x5d1e80c0123206a1ULL,
     0x91b1b3f312ac1be8ULL, 0x19bc3eed09a866a2ULL, 0x29ed711bc0303fb1ULL,
     0x35fa6009ed810bd6ULL, 0x69a90fb1b95d63aeULL, 0x38bedd43ad759c34ULL,
     0xa787e2edc8133e84ULL}, {0x552c119614a1f4e5ULL, 0x68d39af3ec9b831fULL,
     0x3996ac49f898e965ULL, 0x8d79f80825eb64f1ULL, 0x5babe19fdfe4ec70ULL, 
     0x7aa6a18653c9a0ccULL, 0x510e3b69d44d759aULL, 0x3486650c8e18c529ULL},
    {0x2b2b4051acc6e9a5ULL, 0x49a9817e3d56ff74ULL, 0xc6fd82e30671aeb5ULL,
     0x20be895e589409fdULL, 0xd8adec17910a1d4fULL, 0x18882870301f5a0bULL,
     0x225cb56188423cf3ULL, 0x4b61cc25258ea5c8ULL}, {0x63f68eef219e6c01ULL,
     0x7c29a354966d3f55ULL, 0x1cf4c5f3b98ce4deULL, 0x1476b731259de5f8ULL, 
     0x73642b2ec6c2a44bULL, 0x407adff58f6e054cULL, 0x811f7cbc6a7e08c7ULL,
     0x277c0323eb636b90ULL}
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
    
    // declare 8 by multiple of 8 matrix
    matrix = new uint8_t*[8];
    for(int c=0;c<8;c++) matrix[c] = new uint8_t[m_ch];
    
    // add delimeter to denote end of data
    input+=0x90;
    
    // padding
    for(int c=0;c<pad;c++) {
        input+='0';
    }
    
    // add length of input in hex to end of input input
    ss << std::setfill('0') << std::setw(16) << std::hex << len;
    input+=ss.str();
    
    // 1-d array to 2-d matrix
    for(int r=0;r<8;r++) {
        for(int c=0;c<m_ch;c++) {
            matrix[r][c] = input[c*8+r];
        }
    }
    
    return matrix;
}

uint64_t ***prec_kib512(uint8_t **matrix, uint64_t m_ch=8)
{
    // initialize manipulation matrix with input matrix
    uint64_t ***manip_m = nullptr;
    manip_m = new uint64_t **[m_ch/8];
    for(int32_t i=0;i<m_ch/8;i++) {
        manip_m[i] = new uint64_t *[8];
        for(int j=0;j<8;j++) {
            manip_m[i][j] = new uint64_t[8];
            for(int k=0,n=56;k<8,n>=0;k++,n-=8) {
                uint64_t temp=0;
                for(int32_t x=0;x<8;x++) {
                    temp<<=x*2;
                    temp|=matrix[j][x*i];
                }
                std::cout << std::hex << temp << "\n";
                manip_m[i][j][k] = temp << n & 0xffffffffffffffffULL;
                continue;
            }
        }
    }
    // what if there is multiple blocks? has to loop like x=m_ch-8 x < m_ch
    
    if constexpr(std::endian::native == std::endian::little) {
        
    }
    
    return manip_m;
}

int main()
{
    std::string in = "abcdefgh";
    uint8_t **m = kib512_prep(in);
    uint64_t ***manipm = prec_kib512(m);
    
    for(int i=0;i<8;i++) {
        for(int j=0;j<8;j++) {
            // std::cout << std::hex << manipm[0][i][j] << " ";
        }
    }
    
    std::cout << "\n\n" << in.length();
    for(int r=0;r<8;r++) {
        for(int c=0;c<8;c++) {
            // std::cout << m[r][c];
        }
    }
    return 0;
}
