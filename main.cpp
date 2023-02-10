// testing kib512.cpp
#include "kib512.cpp"
#include <thread>
#include <random>
#include <algorithm>
#include <string>
#include <array>
#include <vector>

// benchmark test for proper-diffusion
void test()
{
    // Check for proper diffusion
    std::array<GaloisFieldP*,UINT16_MAX> hashes;
    for(uint32_t i=0;i<UINT16_MAX;i++) {
       // generate random hash
	   std::string str("0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz!@#$%^&*()-_+=[]{};:\'\"/?.,><`~");
	   std::random_device rd; // obtain a random number from hardware
	   std::mt19937 gen(rd()); // seed the generator
	   std::uniform_int_distribution<> distr(0, str.length()); // define the range
	   std::mt19937 generator(rd());
	   std::shuffle(str.begin(), str.end(), generator);
	   std::string rand_str = str.substr(0, distr(gen));
        Kib512 kib512_hash = Kib512(rand_str);
        GaloisFieldP *hash = kib512_hash.hash;
        hashes[i] = hash;
		std::cout << std::endl << "hash: ";
		for(uint8_t j=0;j<8;j++) std::cout << hashes[i][j];

        // check for proper diffusion
        for(uint32_t j=0;j<i;j++) {
			uint8_t diff=0;
			for(uint8_t k=0;k<8;k++) {
				if(hashes[i][k] == hashes[j][k]) diff++;
			}
            if(diff==4) { // if first half is same. It doesn't seem to have proper diffusion
                std::cout << std::endl << "No proper diffusion";
                exit(EXIT_FAILURE);
            }
        }
    }
}

int main(int argc, char **argv) {
    std::string in = "";
    if (argc == 1) {
        in = "abc";
        Kib512 kib512_hash = Kib512(in);
        GaloisFieldP *hash = kib512_hash.hash;
        std::cout << "\n\nhash " << in << ":";
        for(int i=0;i<8;i++) {
            std::cout << std::hex << hash[i] << " ";
        }
    } else {
        for(int i=1;i<argc;i++) {
            in += std::string(argv[i]);
        }
        Kib512 kib512_hash = Kib512(in);
        GaloisFieldP *hash = kib512_hash.hash;
		std::cout << std::endl << "hash of \"" << in << "\": ";
		for(uint8_t i=0;i<8;i++) std::cout << hash[i] << " ";
	}
	test();
	exit(EXIT_FAILURE);
	std::vector<std::thread> threads;
	for(uint8_t i=0;i<10;i++) {
		threads.push_back(std::thread(test));
	}
	for(uint8_t i=0;i<10;i++) {
		threads[i].join();
	}

}
