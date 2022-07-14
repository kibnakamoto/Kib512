#include <iostream>
#include <stdint.h>
#include <cmath>

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
    // try: 11th prime number(37),
    // first prime number (11-7 as index), value=3,
    // 18th prime number 11+7 as index
    for (int i=0;i<64;i+=2) {
        for(int j=0;j<64;j+=3) {
            for(int j=0;j<64;j+=3) {
                
            }
        }
    }
    return 0;
}
