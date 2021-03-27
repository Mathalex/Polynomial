#include <algorithm>
#include <iostream>
#include <map>

using namespace std;

template <typename T>
class Polynomial {
private:
    map<size_t, T> pol;

    void Cut() {
        while (size() > 0 && back() == T(0)) {
            pol.erase(--pol.end());
        }
    }

    void Sieve() {
        for (auto it = pol.begin(); it != pol.end();) {
            if (it->second == T(0)) {
                it = pol.erase(it);
            } else {
                ++it;
            }
        }
    }

    Polynomial& Norm() {
        if (size() > 0) {
            T div = back();
            for (auto& el : pol) {
                el.second /= div;
            }
        }
        return *this;
    }

public:
    Polynomial(const vector<T>& p) {
        for (size_t i = 0; i < p.size(); ++i) {
            if (p[i] != T(0)) {
                pol[i] = p[i];
            }
        }
    }

    Polynomial(const T& p_0 = T()) {
        if (p_0 != T(0)) {
            pol[0] = p_0;
        }
    }

    template <typename Iter>
    Polynomial(Iter first, const Iter& last) {
        for (size_t i = 0; first != last; ++first, ++i) {
            if (*first != T(0)) {
                pol[i] = *first;
            }
        }
    }

    friend bool operator<=>(const Polynomial& p1, const Polynomial& p2) = default;

    friend Polynomial operator+(const Polynomial& p1, const Polynomial& p2) {
        Polynomial res;
        for (const auto& [deg, coef] : p1) {
            res.pol[deg] += coef;
        }
        for (const auto& [deg, coef] : p2) {
            res.pol[deg] += coef;
        }
        res.Sieve();
        return res;
    }

    Polynomial& operator+=(const Polynomial& other) {
        return (*this = *this + other);
    }

    friend Polynomial operator-(const Polynomial& p1, const Polynomial& p2) {
        Polynomial res;
        for (const auto& [deg, coef] : p1) {
            res.pol[deg] += coef;
        }
        for (const auto& [deg, coef] : p2) {
            res.pol[deg] -= coef;
        }
        res.Sieve();
        return res;
    }

    Polynomial& operator-=(const Polynomial& other) {
        return (*this = *this - other);
    }

    friend Polynomial operator*(const Polynomial& p1, const Polynomial& p2) {
        Polynomial res;
        for (const auto& [deg1, coef1] : p1) {
            for (const auto& [deg2, coef2] : p2) {
                res.pol[deg1 + deg2] += coef1 * coef2;
            }
        }
        res.Sieve();
        return res;
    }

    Polynomial& operator*=(const Polynomial& other) {
        return (*this = *this * other);
    }

    T operator[](size_t i) const {
        if (pol.contains(i)) {
            return pol.at(i);
        }
        return T(0);
    }

    int Degree() const {
        if (size() == 0) {
            return -1;
        }
        return pol.rbegin()->first;
    }

    size_t size() const {
        return pol.size();
    }

    const T& back() const {
        return pol.rbegin()->second;
    }

    template<typename Tx>
    Tx FastPow(Tx x, size_t n) const {
        Tx res(1);
        while (n > 0) {
            if (n & 1) {
                res *= x;
            }
            n >>= 1;
            x *= x;
        }
        return res;
    }

    T operator()(const T& x) const {
        T res(0);
        T cur(1);
        size_t prev_deg = 0;
        for (const auto& [deg, coef] : pol) {
            cur *= FastPow<T>(x, deg - prev_deg);
            prev_deg = deg;
            res += cur * coef;
        }
        return res;
    }

    class ConstIterator {
    private:
        typename map<size_t, T>::const_iterator it;

    public:
        ConstIterator(const typename map<size_t, T>::const_iterator& it) : it(it) {}

        friend bool operator<=>(const ConstIterator& p1, const ConstIterator& p2) = default;

        const auto& operator*() const {
            return *it;
        }

        decltype(auto) operator->() const {
            return it.operator->();
        }

        ConstIterator& operator++() {
            ++it;
            return *this;
        }

        ConstIterator operator++(int) {
            ConstIterator res = *this;
            ++it;
            return res;
        }
    };

    ConstIterator begin() const {
        return {pol.cbegin()};
    }

    ConstIterator end() const {
        return {pol.cend()};
    }

    friend ostream& operator<<(ostream& out, const Polynomial& p) {
        if (p.size() == 0) {
            return (out << T(0));
        }
        for (auto it = p.pol.rbegin(); it != p.pol.rend(); ++it) {
            const auto& [i, cur] = *it;
            if (cur == T(0)) {
                continue;
            }
            if (it != p.pol.rbegin() && cur > T(0)) {
                out << '+';
            }
            if (i == 0) {
                out << cur;
                break;
            }
            if (cur == T(-1)) {
                out << '-';
            } else if (cur != T(1)) {
                out << cur << '*';
            }
            out << 'x';
            if (i > 1) {
                out << '^' << i;
            }
        }
        return out;
    }

    Polynomial operator&(const Polynomial& p) const {
        Polynomial res(0);
        Polynomial cur(1);
        size_t prev_deg = 0;
        for (const auto& [deg, coef] : pol) {
            cur *= FastPow<Polynomial>(p, deg - prev_deg);
            prev_deg = deg;
            res += cur * coef;
        }
        return res;
    }

    friend Polynomial operator/(Polynomial p, const Polynomial& q) {
        int deg = p.Degree() - q.Degree();
        if (deg < 0) {
            return T(0);
        }
        Polynomial res;
        while (deg >= 0) {
            T div = res.pol[deg] = p.back() / q.back();
            for (const auto& [i, coef] : q) {
                p.pol[i + deg] -= coef * div;
            }
            p.Cut();
            deg = p.Degree() - q.Degree();
        }
        return res;
    }

    friend Polynomial operator%(const Polynomial& p, const Polynomial& q) {
        return p - p / q * q;
    }

    friend Polynomial operator,(Polynomial p, Polynomial q) {
        while (q != T(0)) {
            Polynomial&& t = p % q;
            p = move(q);
            q = t;
        }
        return p.Norm();
    }
};

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(nullptr);

    vector<double> z(8);
    z[0] = 1;
    z[4] = 1;
    Polynomial<double> f(z.begin(), z.end()), g({-1, -1.5, 1});
    auto it = f.begin();
    cout << (++it)->first << endl;
    cout << f(10) << endl;
    //f *= 3 * g;
    cout << (Polynomial<double>({1, 0, 1}) & f) << endl;
    for (auto& el : f) {
        cout << el.first << ' ' << el.second << endl;
    }
    //cout << (f, g);*/
    return 0;
}
