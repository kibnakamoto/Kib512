#include <iostream>
#include <stdint.h>
#include <string>
#include <sstream>

#if !defined(UINT8_MAX)
    using uint8_t = unsigned char
#endif

// pre-processing of kib512
uint8_t** kib512_prep(std::string input)
{
    uint8_t** matrix = nullptr;
    
    // matrix width and height
    uint64_t len = input.length();
    uint64_t m_wh = ((((512-((len*8)+1)-64) % 512)-7)/8 + len + 9)/8;
    std::cout << m_wh;
    matrix = new uint8_t*[m_wh];
    for(int c=0;c<8;c++) matrix[c] = new uint8_t[m_wh];
    
    // add delimeter to denote end of data
    input += 0x90;
    
    // do not pad if length is 64
    if(input.length() != 64) {
        for(int c=len+1;c<63;c++) {
            input+='0';
        }
    }
    
    // 1-d array to 2-d matrix
    for(int r=0;r<8;r++) {
        for(int c=0;c<8;c++) {
            matrix[r][c] = input[c*8+r];
        }
    }
    
    return matrix;
}

int main()
{
    uint8_t** m = kib512_prep("abcd");
    
    std::cout << "\n\n" << in.length();
    for(int r=0;r<8;r++) {
        for(int c=0;c<8;c++) {
            // std::cout << m[r][c];
        }
    }
    // second delimeter equals length in hex
    return 0;
}
