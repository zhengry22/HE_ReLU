#include<iostream>
using namespace std;

int64_t mod_inverse(int64_t a, int64_t m)
{
    int64_t m0 = m, t, q;
    int64_t x0 = 0, x1 = 1;
    if (m == 1)
        return 0;
    while (a > 1)
    {
        q = a / m;
        t = m;
        m = a % m;
        a = t;
        t = x0;
        x0 = x1 - q * x0;
        x1 = t;
    }
    if (x1 < 0)
        x1 += m0;
    return x1;
}

int main() {
    int inv = mod_inverse(15, 1024);
    cout << inv;
    return 0;
}