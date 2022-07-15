#include <iostream>
#include <stdint.h>
#include <cmath>
#include <bitset>

inline uint64_t rr(uint64_t x, unsigned int n) {
    return (x >> n)|(x << ((sizeof(x)<<3)-n));
}

inline uint64_t lr(uint64_t x, unsigned int n) {
    return (x << n)|(x >> ((sizeof(x)<<3)-n));
}

inline std::string bin(uint64_t x) {
    return std::bitset<64>(x).to_string();
}

int main() {
    uint16_t primes[17];
    uint16_t num = 3;
    int p_index = 0;
    bool isprime;
    while(p_index != sizeof(primes)/sizeof(primes[0])) {
        isprime = true;
        for(int i=2;i<=num/2;i++) {
            if(num%i == 0) {
                isprime = false;
                break;
            }
        }
        if(isprime) {
            primes[p_index] = num;
            p_index++;
        }
        num++;
    }
    
    for(int i=0;i<17;i++) std::cout << primes[i] << " ";
    std::cout << std::endl << std::endl;
    // try: 11th prime number(37),
    // try: 1st prime number(3),
    // try: 16th prime number(59)
    // try: try 9th prime number(29)
    
    uint64_t i=primes[10];
    uint64_t j=primes[0];
    uint64_t k=primes[15];
    uint64_t n=primes[8];
    
    // test rotation on single value(test num)
    uint64_t tn1 = 0x159a300c24f5dc1eULL;
    uint64_t tn2 = 0x0000000000000000ULL;//0x3030303030303030ULL;
    uint64_t tn3 = 0x0000000000000000ULL;//0x3030303030303030ULL;
    uint64_t tn4 = 0x0000000000000000ULL;
    
    uint64_t finalA = rr(tn1, i) xor lr(tn2, j) xor (tn3 << k) xor (tn4 << n) & 0xffffffffffffffffULL;
    uint64_t finalB = rr(tn4, n) xor lr(tn1, i) xor (tn2 << j) xor (tn3 << k) & 0xffffffffffffffffULL;
    uint64_t finalC = rr(tn3, k) xor lr(tn4, n) xor (tn1 << i) xor (tn2 << j) & 0xffffffffffffffffULL;
    uint64_t finalD = rr(tn2, j) xor lr(tn3, k) xor (tn4 << n) xor (tn1 << i) & 0xffffffffffffffffULL;
    uint64_t final = (finalA + finalB + finalC + finalD) & 0xffffffffffffffffULL;
    
    std::cout << std::hex << bin(final) << std::endl << std::endl;
    std::cout << std::hex << bin(tn1) << std::endl;
    std::cout << std::hex << bin(tn2) << std::endl;
    std::cout << std::hex << bin(tn3) << std::endl;
    return 0;
}
