#include <iostream>
#include <stdint.h>
#include <string>

#if !defined(UINT8_MAX)
    using uint8_t = unsigned char
#endif

// pre-processing of kib512
uint8_t** kib512_prepr(std::string input)
{
    uint8_t** matrix = nullptr;
    matrix = new uint8_t*[8];
    uint64_t len = input.length();
    for(int c=0;c<8;c++) matrix[c] = new uint8_t[8];
    
    // add delimeter to denote end of data
    input[len-1] = ;
    
    // do not pad if length is 64
    if(input.length() != 64) {
        for(int c=len+1;c<63;c++) {
            input+='0';
        }
    }
    
    // 1-d array to 2-d matrix
    for(int r=0;r<8;r++) {
        for(int c=0;c<8;c++) {
            matrix[r][c] = input[r*8+c];
        }
    }
    
    return matrix;
}

int main()
{
    uint8_t** m = kib512_prepr("abcd");
    for(int r=0;r<8;r++) {
        for(int c=0;c<8;c++) {
            std::cout << m[r][c];
        }
    }
    
    return 0;
}
