#include <iostream>
#include <stdint.h>
#include <string>

#if !defined(UINT8_MAX)
    using uint8_t = unsigned char
#endif

// pre-processing of kib512
uint8_t** kib512_prep(std::string input)
{
    uint8_t** matrix = nullptr;
    matrix = new uint8_t*[8];
    uint64_t len = input.length();
    for(int c=0;c<8;c++) matrix[c] = new uint8_t[8];
    
    // add delimeter to denote end of data
    input += 0x90>>4;
    
    // do not pad if length is 64
    if(input.length() != 64) {
        for(int c=len+1;c<63;c++) {
            input+='0'-48;
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
    for(int r=0;r<8;r++) {
        for(int c=0;c<8;c++) {
            // std::cout << m[r][c];
        }
    }
    uint8_t sigmax1d = 0x90 >> 4; // sigmax as 1-d
    return 0;
}
