#ifndef __RANDOM_H__
#define __RANDOM_H__
typedef unsigned long int U_Int_32_t;
typedef unsigned long long int  U_Int_64_t;
class Random {
public:
    U_Int_64_t u,v,w;
    Random(U_Int_64_t j) : v(4101842887655102017LL), w(1) {
        u = j ^ v; int64();
        v = u; int64();
        w = v; int64();
    }
    inline U_Int_64_t int64() {
        u = u * 2862933555777941757LL + 7046029254386353087LL;
        v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
        w = 4294957665U*(w & 0xffffffff) + (w >> 32);
        U_Int_64_t x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
        return (x + v) ^ w;
    }
    inline double doub() { return 5.42101086242752217E-20 * int64(); }
    inline unsigned int int32() { return (unsigned int)int64(); }
};
#endif
