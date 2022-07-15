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

/*
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
         0xa787e2edc8133e84ULL}
*/

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
    uint64_t n=1;//primes[8];
    
    // test rotation on single value(test num)
    uint64_t tn1 = 0x159a300c24f5dc1eULL;
    uint64_t tn2 = 0x3030303030303030ULL;
    uint64_t tn3 = 0x3030303030303030ULL;
    uint64_t tn4 = 0x0000000000000000ULL;
    
    uint64_t finalA = rr(tn1, i) xor lr(tn2, j) xor (tn3 << k) xor (tn4 << n) & 0xffffffffffffffffULL;
    uint64_t finalB = rr(tn4, n) xor lr(tn1, i) xor (tn2 << j) xor (tn3 << k) & 0xffffffffffffffffULL;
    uint64_t finalC = rr(tn3, k) xor lr(tn4, n) xor (tn1 << i) xor (tn2 >> j) & 0xffffffffffffffffULL;
    uint64_t finalD = rr(tn2, j) xor lr(tn3, k) xor (tn4 << n) xor (tn1 >> i) & 0xffffffffffffffffULL;
    uint64_t final = (finalA + finalB + finalC + finalD) & 0xffffffffffffffffULL;
    
    std::cout << std::hex << final << std::endl;
    std::cout << std::hex << bin(finalA) << std::endl;
    std::cout << std::hex << bin(finalB) << std::endl;
    std::cout << std::hex << bin(finalC) << std::endl;
    std::cout << std::hex << bin(finalD) << std::endl << std::endl;

    std::cout << std::hex << bin(tn1) << std::endl;
    std::cout << std::hex << bin(tn2) << std::endl;
    std::cout << std::hex << bin(tn3) << std::endl;
    return 0;
}
