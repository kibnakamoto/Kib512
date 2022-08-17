#include <iostream>
#include <stdint.h>
#include <cmath>
#include <bitset>
#include <string>
#include <sstream>
#include <iomanip>
#include <bit>
#include <assert.h>

#if !defined(UINT8_MAX)
using uint8_t = unsigned char
#endif

#if !defined(UINT32_MAX)
using uint32_t = unsigned int
#endif

#if !defined(UINT64_MAX)
using uint64_t = unsigned long long
#endif

template<size_t size> class Matrix;

// Operators in Galois Field
class GaloisFieldP {
    public:
    uint64_t x;
    uint64_t p;
    
    GaloisFieldP(uint64_t integer=0, uint64_t prime=0xffffffffffffffc5ULL) {
        if (prime == 0) {
            throw std::invalid_argument("p = 0, cannot do modulo 0");
        }
        this->x = integer;
        this->p = prime;
    }
    
    // modular inverse
    GaloisFieldP operator~() const;
    
    GaloisFieldP operator= (GaloisFieldP const &y) {
        p = y.p;
        x = y.x%p;
        return *this;
    }
    
    // assignment operator for integers
    GaloisFieldP operator= (uint64_t &y) {
        x = y%p;
        return *this;
    }
    
    GaloisFieldP operator+ (GaloisFieldP const &y) {
        GaloisFieldP s(((__uint128_t)x + y.x)%p, p);
        return s;
    }
    
    GaloisFieldP operator+ (const uint64_t &y) {
        GaloisFieldP s(((__uint128_t)x + y)%p, p);
        return s;
    }
    
    GaloisFieldP operator- (const uint64_t &y) {
        GaloisFieldP s(((__uint128_t)x - y+p)%p, p);
        return s;
    }
    
    GaloisFieldP operator- (GaloisFieldP const &y) {
        GaloisFieldP s(((__uint128_t)x - y.x+p)%p, p);
        return s;
    }
    
    GaloisFieldP operator-= (GaloisFieldP const &y) {
        x = ((__uint128_t)x - y.x+p)%p;
        return *this;
    }
    
    GaloisFieldP operator+= (GaloisFieldP const &y) {
        x = ((__uint128_t)x + y.x)%p;
        return *this;
    }
    
    GaloisFieldP operator+= (const uint64_t &y) {
        x = ((__uint128_t)x + y)%p;
        return *this;
    }
    
    GaloisFieldP operator% (GaloisFieldP const &y) {
        if (y.x != 0) {
            GaloisFieldP s(x%(y.x), p);
            return s;
        }
        return 0;
    }
    
    // modulo for integer type, for overloading modulo for % y.p instead of % y.x
    GaloisFieldP operator% (const uint64_t &y) {
        GaloisFieldP s(x%y, p);
        return s;
    }
    
    GaloisFieldP operator%= (const uint64_t &y) {
        x%= y;
        return *this;
    }
    
    GaloisFieldP operator%= (GaloisFieldP const &y) {
        x%= y.x;
        return *this;
    }
    
    GaloisFieldP operator<<= (GaloisFieldP const &y) {
        x = (x << y.x)%p;
        return *this;
    }
    
    GaloisFieldP operator<< (GaloisFieldP const &y) {
        GaloisFieldP s((x << y.x)%p);
        return s;
    }
    
    GaloisFieldP operator<< (uint64_t &y) {
        GaloisFieldP s((x << y)%p, p);
        return s;
    }

    
    GaloisFieldP operator>>= (GaloisFieldP const &y) {
        x = (x >> y.x)%p;
        return *this;
    }
    
    GaloisFieldP operator>> (GaloisFieldP const &y) {
        GaloisFieldP s((x >> y.x)%p, p);
        return s;
    }
    
    GaloisFieldP operator>> (uint64_t const &y) {
        GaloisFieldP s((x >> y)%p, p);
        return s;
    }
    
    GaloisFieldP operator| (uint64_t const &y) {
        GaloisFieldP s((x | y)%p, p);
        return s;
    }
    
    GaloisFieldP operator| (GaloisFieldP const &y) {
        GaloisFieldP s((x | y.x)%p, p);
        return s;
    }
    
    GaloisFieldP operator|= (GaloisFieldP const &y) {
        x = (x | y.x)%p;
        return *this;
    }
    
    bool operator> (GaloisFieldP const &y) {
        bool s = x > y.x;
        return s;
    }
    
    bool operator< (GaloisFieldP const &y) {
        bool s = x < y.x;
        return s;
    }
    
    bool operator== (GaloisFieldP const &y) {
        bool s = x == y.x;
        return s;
    }
    
    GaloisFieldP operator& (GaloisFieldP const &y) {
        GaloisFieldP s(x & y.x, p);
        return s;
    }
    
    bool operator&& (GaloisFieldP const &y) {
        bool s = x & y.x;
        return s;
    }

    
    GaloisFieldP operator&= (GaloisFieldP const &y) {
        x = (x & y.x)%p;
        return *this;
    }
    
    GaloisFieldP operator^ (GaloisFieldP const &y) {
        GaloisFieldP s((x ^ y.x)%p, p);
        return s;
    }
    
    GaloisFieldP operator^= (GaloisFieldP const &y) {
        x = (x ^ y.x)%p;
        return *this;
    }
    
    GaloisFieldP operator/ (GaloisFieldP const &y) {
        GaloisFieldP s((x / y.x)%p, p);
        return s;
    }
    
    GaloisFieldP operator/= (GaloisFieldP const &y) {
        x = (x / y.x)%p;
        return *this;
    }
    
    GaloisFieldP operator* (const uint64_t &y) {
        GaloisFieldP mul = GaloisFieldP(((__uint128_t)x*y)%p,p);
        return mul;
    }
    
    GaloisFieldP operator* (GaloisFieldP y) {
        GaloisFieldP mul = GaloisFieldP(((__uint128_t)x*y.x)%p,p);
        return mul;
    }
    
    GaloisFieldP operator*= (GaloisFieldP y) {
        x = ((__uint128_t)x*y.x)%p;
        return *this;
    }
    
    GaloisFieldP operator*= (uint64_t y) {
        x = ((__uint128_t)x*y)%p;
        return *this;
    }
};

// output stream operator
std::ostream& operator<< (std::ostream& out, GaloisFieldP toprint) {
    out << toprint.x;
    return out;
}

using point_t = std::pair<GaloisFieldP,GaloisFieldP>;

// extended Eucludian Algorithm
inline GaloisFieldP extended_gcd(GaloisFieldP a, GaloisFieldP b, GaloisFieldP &x)
{
    x = 1;
    GaloisFieldP y = GaloisFieldP(0,b.p);
    if (b == 0) return a;
    GaloisFieldP n_x = GaloisFieldP(0,b.p);
    GaloisFieldP n_y = GaloisFieldP(1,b.p);
    GaloisFieldP n_r = GaloisFieldP(b.p,0xffffffffffffffffULL);
    GaloisFieldP r = GaloisFieldP(a.x,b.p); 
    GaloisFieldP quot, tmp;
    while ((bool)n_r.x) {
        quot = r / n_r;
        tmp = r;
        r = n_r;
        n_r = tmp - quot * n_r;
        tmp = x;
        x = n_x;
        n_x = tmp - quot * n_x;
        tmp = y;
        y = n_y;
        n_y = tmp - quot * n_y;
    }
    return r;
}

inline GaloisFieldP mod_inv(GaloisFieldP a, GaloisFieldP p) {
    GaloisFieldP x;
    GaloisFieldP gcd = extended_gcd(a,p,x);
    if (gcd.x != 1) {
        std::cout << std::flush << "gcd isn't one\n";
        return 1;
    }
    return (x+p) % p;
}

// modular inverse
GaloisFieldP GaloisFieldP::operator~() const {
    GaloisFieldP s(mod_inv(x,p).x, p);
    return s;
}

// only for square matrix multipication on GF(p) for same size matrices
template<size_t size>
class Matrix
{
    public:
    GaloisFieldP **m1;
    GaloisFieldP p;
    GaloisFieldP **res;
    
    Matrix(GaloisFieldP matrix[size][size], GaloisFieldP prime) {
        m1 = nullptr;
        m1 = new GaloisFieldP*[size];
        for(size_t i=0;i<size;i++) {
            m1[i] = new GaloisFieldP[size];
            
            // copy parameter matrix to operate on
            for(size_t j=0;j<size;j++) {
                m1[i][j] = matrix[i][j] % prime;
            }
        }
        p = prime;
    }
    
    // square matrix multipication for 2 2d matrices on GF(p)
    Matrix operator* (GaloisFieldP m[size][size]) {
        res = nullptr;
        res = new GaloisFieldP*[size];
        for(size_t i=0;i<size;i++) {
            res[i] = new GaloisFieldP[size];
            for(size_t j=0;j<size;j++) {
                res[i][j].p = p.p;  // change field sizes
                for(size_t k=0;k<size;k++) {
                    m[i][k].p = p.p; // change field sizes
                    res[i][j] = res[i][j] + m[i][k] * m1[k][j];
                }
            }
        }
        return *this;
    }
    
    GaloisFieldP* operator[] (int64_t index) {
        return res[index];
    }
    
    // ~Matrix() {
        // for(int i=0;i<size;i++) {
            // delete[] res[i];
            // delete[] m1[i];
        // }
        // delete[] res;
        // delete[] m1;
    // }
};

// Taha Canturk Kibnakamoto 64-bit Koblitz curve with prime field size
struct Tckp64k1
{
    const GaloisFieldP a = 0x0000000000000000ULL; // taken from SEC curves domain parameters
    const GaloisFieldP b = 0x0000000000000007ULL; // taken from SEC curves domain parameters
    
    // generated without the use of SEC specifications on how to generate p and q
    // since that is generated randomly, there isn't a need to
    // verified as prime using Fermat's Little Theorem in GF(gf_p) where gf_p
    // denotes the potential largest unsigned 64-bit prime number
    const GaloisFieldP p = 0xffffffffffffffc5ULL; // largest unsigned 64-bit prime number
    
    // calculated with sagemath
    // __uint128_t n = 0x100000001dc431000; // bigger than field size
    
    // calculated with sagemath
    const GaloisFieldP gx = 0x7143332d09966ea9ULL; // x coordinate of generator point
    const GaloisFieldP gy = 0xb5fbd04e6e22f933ULL; // y coordinate of generator point
    const GaloisFieldP h = 0x0000000000000001ULL; // co-factor
};

// bitwise right-rotate
inline GaloisFieldP rr(GaloisFieldP x, unsigned int n) { return (x >> n)|(x << (64-n)); }

// bitwise left-rotate
inline GaloisFieldP lr(GaloisFieldP x, unsigned int n) { return (x << n)|(x >> (64-n)); }

// convert to bit string
inline std::string bin(uint64_t x) { return std::bitset<64>(x).to_string(); }

inline point_t point_add(point_t p1, point_t p2, GaloisFieldP p, GaloisFieldP a)
{
    // equation for calculating point addition in ECC
    // find lambda
    
    GaloisFieldP xp,yp,xq,yq,__lambda;
    xp = std::get<0>(p1);
    yp = std::get<1>(p1);
    xq = std::get<0>(p2);
    yq = std::get<1>(p2);
    if (yp == yq || xp == xq) {
        __lambda = ((GaloisFieldP(3,p.p)*xp*xp + a)*~(GaloisFieldP(2,p.p)*yp)) % p;
    } else {
        __lambda = (yq-yp)*~(xq-xp);
    }
    GaloisFieldP xr = (__lambda*__lambda - xp - xq)%p;
    GaloisFieldP yr = (__lambda*(xp-xr) - yp) % p;
    
    /* C++ modulo works differently than python's so subtract 1 from final output */
    
    return std::make_pair(xr+p, yr);
}

inline point_t point_double(point_t p1, GaloisFieldP p, GaloisFieldP a) {
    GaloisFieldP x,y;
    x = std::get<0>(p1);
    y = std::get<1>(p1);
    GaloisFieldP __lambda = (GaloisFieldP(3,p.p)*(x*x) + a)*~(GaloisFieldP(2,p.p)*y);
    std::cout << ~(GaloisFieldP(2,p.p)*y);
    GaloisFieldP xr = __lambda*__lambda - GaloisFieldP(2,p.p)*x;
    GaloisFieldP yr = __lambda*(x - xr) - y;
    return std::make_pair(xr, yr);
}

// montgomery ladder ECC point multiplication 
inline point_t montgomery_ladder(point_t r0, uint64_t k, uint64_t p,
                                 uint64_t a) {
    point_t r1 = point_double(r0,p,a);
    std::string bits = bin(k).erase(0);
    for(uint8_t i : bits) {
        if (i == '0') {
            r1 = point_add(r0,r1,p,a);
            r0 = point_double(r0,p,a);
        } else {
            r0 = point_add(r0,r1,p,a);
            r1 = point_double(r1,p,a);
        }
    }
    return r0;
}

// pre-processing of kib512
class Kib512 {
    public: // TODO: make private after debugging
        Tckp64k1 curve; // curve domain parameters
        uint8_t **matrix;
        uint64_t ***manip_m;
        
        // amount of blocks to be processed
        uint64_t b_size;
    
    public:
    uint64_t m_ch;
    static constexpr uint64_t const const_m[8][8] = {
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
    
    mutable uint64_t hash[8] {
        0x5287173768046659ULL, 0x1388c1a81885db29ULL, 0xc582055c7b0f1a24ULL,
        0x9be22d058c2ae082ULL, 0xa36b2344c3b2e0d0ULL, 0x8b74c2c08074d4b1ULL,
        0x2a1c87eb1011bd80ULL, 0x64edcba54a4d0a5aULL
    };
    
    // galois field size
    GaloisFieldP gf_p = curve.p; // prime field size
    
    Kib512(std::string input);
    
    ~Kib512();
    
    void kib512_prep(std::string input)
    {
        matrix = nullptr;
        std::stringstream ss;
        
        // matrix column height
        uint64_t len = input.length();
        uint32_t pad = (((512-((len*8)+1)-128) % 512)-7)/8;
        m_ch = (pad + len + 17)/8;
        
        // declare 8 by multiple of 8 matrix
        matrix = new uint8_t*[8];
        for(int c=0;c<8;c++) matrix[c] = new uint8_t[m_ch];
        
        // add delimeter to denote end of data
        input+=0x90;
        
        // padding
        for(int c=0;c<pad;c++) input+='0';
        
        // add length of input in hex to end of input input
        ss << std::setfill('0') << std::setw(16) << std::hex << len;
        input+=ss.str();
        
        // 1-d array to 2-d matrix
        for(int r=0;r<8;r++) {
            for(uint64_t c=0;c<m_ch;c++) {
                matrix[r][c] = (uint8_t)input[r*m_ch+c]; // input ordered as row
            }
        }
        
        // for(int r=0;r<8;r++) {
        //     for(uint64_t c=0;c<m_ch;c++) {
        //         std::cout << std::hex << matrix[r][c]+0;
        //     }
        // }std::cout << std::endl;
    }
    
    void prec_kib512()
    {
        b_size = m_ch/8;
        
        // initialize manipulation matrix with input matrix
        manip_m = nullptr;
        manip_m = new uint64_t **[b_size];
        uint64_t loop_count = 0;
        
        // initialize 3-d matrix to zero to indexes that won't be initialized to
        // values of 2-d matrix
        for(int i=0;i<b_size;i++) {
            manip_m[i] = new uint64_t *[8];
            manip_m[i][0] = new uint64_t[8];
            for(int j=1;j<8;j++) { // to avoid unnecesary padding, j=1
                manip_m[i][j] = new uint64_t[8];
                for(int k=0;k<8;k++) manip_m[i][j][k] = 0x0000000000000000ULL;
            }
        }
        
        uint64_t mmi=0; // manipulation matrix i
        int mmj=0;  // manipulation matrix j
        for(int j=0;j<8;j++) {
            for(uint64_t i=0;i<b_size;i++) {
                GaloisFieldP temp=0;
                // add matrix values while avoiding repetition of values
                for(int x=0;x<8;x++) {
                    temp or_eq (GaloisFieldP)matrix[j][x+i*8] << 56-x*8;
                }
                
                    manip_m[mmi][0][mmj] = temp.x;
                mmj = (mmj+1)%8;
            }
            if(mmj%8==0) mmi++;
        }
        
        // pre-compression. Get rid of extra padding
        for(uint64_t i=0;i<b_size;i++) {
            for(int j=1;j<8;j++) {
                for(int k=0;k<8;k++) {
                    // primes used for rotation and shifting
                    // chosen so that one big, one small one rotation is made each time
                    unsigned int p[4] = {37, 3, 59, 5};
                    
                    // matrix indexes are chosen within reason, +7,+6 is for
                    // length of matrix and +1 and +0 is for start of the message.
                    
                    GaloisFieldP tn1,tn2,tn3,tn4 = manip_m[(i+b_size)%b_size][j-1][(k+7)%8];
                    tn3 = manip_m[(i+b_size)%b_size][j-1][(k+6)%8];
                    tn2 = manip_m[i][j-1][(k+1)%8];
                    tn1 = manip_m[i][j-1][k];
                    GaloisFieldP sigma0 = rr(tn1, p[0]) xor rr(tn2, p[1]) xor (tn3 << p[2]) |
                                      (tn4 << p[3]) % gf_p;
                    GaloisFieldP sigma1 = rr(tn4, p[3]) xor lr(tn1, p[0]) xor (tn2 << p[1]) |
                                      (tn3 << p[2]) % gf_p;
                    GaloisFieldP sigma2 = rr(tn3, p[2]) xor lr(tn4, p[3]) xor (tn1 << p[0]) |
                                      (tn2 << p[1]) % gf_p;
                    GaloisFieldP sigma3 = rr(tn2, p[1]) xor lr(tn3, p[2]) xor (tn4 << p[3]) |
                                      (tn1 << p[0]) % gf_p;
                    
                    // use arithmetic addition for non-linearity
                    manip_m[i][j][k] = (sigma0/p[k%4] +
                                        sigma1 + sigma2 + sigma3).x;
                }
            }
        }
    }
    
    // compression function
    void hash_kib512()
    {
        
    }
    
    template<typename T>
    T hashstr(std::string input) {
        // check if type of T is valid, then return hash
        assert((typeid(T) != typeid(std::string) || typeid(T) !=
                typeid(uint64_t*)) && "Kib512 hashstr:\t wrong type detected");
        
        // calculate hash
        kib512_prep(input);
        prec_kib512();
        hash_kib512();
        
        // return hash as requested data type
        if(typeid(T) == typeid(std::string)) {
            std::stringstream ss;
            for(int i=0;i<8;i++) {
                ss << std::setfill('0') << std::setw(16) << std::hex << hash[i];
            }
            return ss.str();
        }
        else {
            return hash;
        }
    }
};

// define constructor
Kib512::Kib512(std::string input) {
    
}

// define destructor
Kib512::~Kib512() {
    for(uint64_t i=0;i<b_size;i++) {
        for(int j=0;j<8;j++) {
            delete[] manip_m[i][j];
        }
        delete[] matrix[i];
    }
    delete[] matrix;
    for(int i=0;i<b_size;i++) {
        delete[] manip_m[i];
    }
    delete[] manip_m;
}

int main() {
    Kib512 *kib512 = new Kib512("input");
    // std::string in = "abcdefghqwertyuioplkjhgfdsazxcvbnm1234567890!@#$%^&*()\\/";
    std::string in = "abc";
    kib512->kib512_prep(in);
    kib512->prec_kib512();
    
    for(int i=0;i<kib512->b_size;i++) {
        for(int j=0;j<8;j++) {
            for(int k=0;k<8;k++) {
                std::cout << std::setfill('0') << std::setw(16) << std::hex
                          << kib512->manip_m[i][j][k] << "\t";
            }
        }
    }
    std::cout << "\n\n" << std::hex << in.length();
    return 0;
}
