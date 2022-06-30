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
    
    // get 64-bit hex length in hex
    ss << std::setfill('0') << std::setw(16) << std::hex << (len<<2);
    std::string hexlen = ss.str();
    input+=hexlen;
    std::cout  << std::setfill('0') << std::setw(16) << std::hex << (len*4);
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
    // second delimeter equals length in hex
    return 0;
}
