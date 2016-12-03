// Input procedures
// Copyright (c) 2016 Ralf Quast
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
#ifndef ESPECIA_READLINE_H
#define ESPECIA_READLINE_H

#include <cstddef>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace especia {
    template<class A>
    std::istream &read(std::istream &is, std::vector<A> &a, size_t n, bool append = false);

    template<class A, class B>
    std::istream &read(std::istream &is, std::vector<A> &a, std::vector<B> &b, size_t n, bool append = false);

    template<class A, class B, class C>
    std::istream &read(std::istream &is, std::vector<A> &a, std::vector<B> &b,
                       std::vector<C> &c, size_t n, bool append = false);

    template<class A, class B, class C, class D>
    std::istream &read(std::istream &is, std::vector<A> &a, std::vector<B> &b,
                       std::vector<C> &c, std::vector<D> &d, size_t n, bool append = false);

    template<class A, class B, class C, class D, class E>
    std::istream &read(std::istream &is, std::vector<A> &a, std::vector<B> &b,
                       std::vector<C> &c, std::vector<D> &d, std::vector<E> &e, size_t n,
                       bool append = false);

    template<class A, class B, class C, class D>
    std::istream &read(std::istream &is, std::vector<A> &a, std::vector<B> &b,
                       std::vector<C> &c, std::vector<D> &d, std::vector<std::string> &s, size_t n,
                       char eol, bool append = false);

    // Read a line, stripping empty lines and comments
    std::istream &readline(std::istream &is, std::string &s, char comment_mark = 0, char eol = '\n',
                           bool eat_empty = true);
}

template<class A>
std::istream &especia::read(std::istream &is, std::vector<A> &a, size_t n, bool append) {
    using namespace std;

    vector<A> ta;

    ta.reserve(n);

    for (size_t i = 0; i < n; ++i) {
        A aa;

        if (is >> aa)
            ta.push_back(aa);
    }

    if (is) {
        if (append) {
            a.insert(a.end(), ta.begin(), ta.end());
        } else {
            a.assign(ta.begin(), ta.end());
        }
    }
    return is;
}

template<class A, class B>
std::istream &especia::read(std::istream &is, std::vector<A> &a, std::vector<B> &b, size_t n, bool append) {
    using namespace std;

    vector<A> ta;
    vector<B> tb;

    ta.reserve(n);
    tb.reserve(n);

    for (size_t i = 0; i < n; ++i) {
        A aa;
        B bb;

        if (is >> aa >> bb) {
            ta.push_back(aa);
            tb.push_back(bb);
        }
    }

    if (is) {
        if (append) {
            a.insert(a.end(), ta.begin(), ta.end());
            b.insert(b.end(), tb.begin(), tb.end());
        } else {
            a.assign(ta.begin(), ta.end());
            b.assign(tb.begin(), tb.end());
        }
    }
    return is;
}

template<class A, class B, class C>
std::istream &especia::read(std::istream &is, std::vector<A> &a, std::vector<B> &b, std::vector<C> &c, size_t n, bool append) {
    using namespace std;

    vector<A> ta;
    vector<B> tb;
    vector<C> tc;

    ta.reserve(n);
    tb.reserve(n);
    tc.reserve(n);

    for (size_t i = 0; i < n; ++i) {
        A aa;
        B bb;
        C cc;

        if (is >> aa >> bb >> cc) {
            ta.push_back(aa);
            tb.push_back(bb);
            tc.push_back(cc);
        }
    }

    if (is) {
        if (append) {
            a.insert(a.end(), ta.begin(), ta.end());
            b.insert(b.end(), tb.begin(), tb.end());
            c.insert(c.end(), tc.begin(), tc.end());
        } else {
            a.assign(ta.begin(), ta.end());
            b.assign(tb.begin(), tb.end());
            c.assign(tc.begin(), tc.end());
        }
    }
    return is;
}

template<class A, class B, class C, class D>
std::istream &especia::read(std::istream &is, std::vector<A> &a, std::vector<B> &b, std::vector<C> &c,
         std::vector<D> &d, size_t n, bool append) {
    using namespace std;

    vector<A> ta;
    vector<B> tb;
    vector<C> tc;
    vector<D> td;

    ta.reserve(n);
    tb.reserve(n);
    tc.reserve(n);
    td.reserve(n);

    for (size_t i = 0; i < n; ++i) {
        A aa;
        B bb;
        C cc;
        D dd;

        if (is >> aa >> bb >> cc >> dd) {
            ta.push_back(aa);
            tb.push_back(bb);
            tc.push_back(cc);
            td.push_back(dd);
        }
    }

    if (is) {
        if (append) {
            a.insert(a.end(), ta.begin(), ta.end());
            b.insert(b.end(), tb.begin(), tb.end());
            c.insert(c.end(), tc.begin(), tc.end());
            d.insert(d.end(), td.begin(), td.end());
        } else {
            a.assign(ta.begin(), ta.end());
            b.assign(tb.begin(), tb.end());
            c.assign(tc.begin(), tc.end());
            d.assign(td.begin(), td.end());
        }
    }
    return is;
}

template<class A, class B, class C, class D, class E>
std::istream &especia::read(std::istream &is, std::vector<A> &a, std::vector<B> &b, std::vector<C> &c, std::vector<D> &d,
         std::vector<E> &e, size_t n, bool append) {
    using namespace std;

    vector<A> ta;
    vector<B> tb;
    vector<C> tc;
    vector<D> td;
    vector<E> te;

    ta.reserve(n);
    tb.reserve(n);
    tc.reserve(n);
    td.reserve(n);
    te.reserve(n);

    for (size_t i = 0; i < n; ++i) {
        A aa;
        B bb;
        C cc;
        D dd;
        E ee;

        if (is >> aa >> bb >> cc >> dd >> ee) {
            ta.push_back(aa);
            tb.push_back(bb);
            tc.push_back(cc);
            td.push_back(dd);
            te.push_back(ee);
        }
    }

    if (is) {
        if (append) {
            a.insert(a.end(), ta.begin(), ta.end());
            b.insert(b.end(), tb.begin(), tb.end());
            c.insert(c.end(), tc.begin(), tc.end());
            d.insert(d.end(), td.begin(), td.end());
            e.insert(e.end(), te.begin(), te.end());
        } else {
            a.assign(ta.begin(), ta.end());
            b.assign(tb.begin(), tb.end());
            c.assign(tc.begin(), tc.end());
            d.assign(td.begin(), td.end());
            e.assign(te.begin(), te.end());
        }
    }
    return is;
}

template<class A, class B, class C, class D>
std::istream &especia::read(std::istream &is, std::vector<A> &a, std::vector<B> &b, std::vector<C> &c, std::vector<D> &d,
         std::vector<std::string> &s, size_t n,
         char eol, bool append) {
    using namespace std;

    vector<A> ta;
    vector<B> tb;
    vector<C> tc;
    vector<D> td;
    vector<string> ts;

    ta.reserve(n);
    tb.reserve(n);
    tc.reserve(n);
    td.reserve(n);
    ts.reserve(n);

    for (size_t i = 0; i < n; ++i) {
        A aa;
        B bb;
        C cc;
        D dd;
        string ss;

        if (is >> aa >> bb >> cc >> dd and getline(is, ss, eol)) {
            istringstream ist(ss);

            ist >> ss;
            if (!ist)
                ss.erase();

            ta.push_back(aa);
            tb.push_back(bb);
            tc.push_back(cc);
            td.push_back(dd);
            ts.push_back(ss);
        }
    }

    if (is) {
        if (append) {
            a.insert(a.end(), ta.begin(), ta.end());
            b.insert(b.end(), tb.begin(), tb.end());
            c.insert(c.end(), tc.begin(), tc.end());
            d.insert(d.end(), td.begin(), td.end());
            s.insert(s.end(), ts.begin(), ts.end());
        } else {
            a.assign(ta.begin(), ta.end());
            b.assign(tb.begin(), tb.end());
            c.assign(tc.begin(), tc.end());
            d.assign(td.begin(), td.end());
            s.assign(ts.begin(), ts.end());
        }
    }
    return is;
}

#endif // ESPECIA_READLINE_H
