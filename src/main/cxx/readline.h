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
#ifndef RQ_READLINE_H
#define RQ_READLINE_H

#include <cstddef>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace RQ {
    template<class A>
    std::istream &read(std::istream &is, std::vector<A> &a, size_t n,
                       bool append = false);

    template<class A, class B>
    std::istream &read(std::istream &is, std::vector<A> &a, std::vector<B> &b, size_t n,
                       bool append = false);

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

    // Read a line while skipping empty lines and stripping comments
    std::istream &readline(std::istream &is, std::string &s,
                           char comment_mark = 0, char eol = '\n');
}

template<class A>
std::istream &
RQ::read(std::istream &is, std::vector<A> &a, size_t n, bool append) {
    using namespace std;

    vector<A> ta;

    ta.reserve(n);

    for (size_t i = 0; i < n; ++i) {
        A a;

        if (is >> a)
            ta.push_back(a);
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
std::istream &
RQ::read(std::istream &is, std::vector<A> &a, std::vector<B> &b, size_t n,
         bool append) {
    using namespace std;

    vector<A> ta;
    vector<B> tb;

    ta.reserve(n);
    tb.reserve(n);

    for (size_t i = 0; i < n; ++i) {
        A a;
        B b;

        if (is >> a >> b) {
            ta.push_back(a);
            tb.push_back(b);
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
std::istream &
RQ::read(std::istream &is, std::vector<A> &a, std::vector<B> &b, std::vector<C> &c,
         size_t n, bool append) {
    using namespace std;

    vector<A> ta;
    vector<B> tb;
    vector<C> tc;

    ta.reserve(n);
    tb.reserve(n);
    tc.reserve(n);

    for (size_t i = 0; i < n; ++i) {
        A a;
        B b;
        C c;

        if (is >> a >> b >> c) {
            ta.push_back(a);
            tb.push_back(b);
            tc.push_back(c);
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
std::istream &
RQ::read(std::istream &is, std::vector<A> &a, std::vector<B> &b, std::vector<C> &c,
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
        A a;
        B b;
        C c;
        D d;

        if (is >> a >> b >> c >> d) {
            ta.push_back(a);
            tb.push_back(b);
            tc.push_back(c);
            td.push_back(d);
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
std::istream &
RQ::read(std::istream &is, std::vector<A> &a, std::vector<B> &b, std::vector<C> &c,
         std::vector<D> &d, std::vector<E> &e, size_t n, bool append) {
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
        A a;
        B b;
        C c;
        D d;
        E e;

        if (is >> a >> b >> c >> d >> e) {
            ta.push_back(a);
            tb.push_back(b);
            tc.push_back(c);
            td.push_back(d);
            te.push_back(e);
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
std::istream &
RQ::read(std::istream &is, std::vector<A> &a, std::vector<B> &b, std::vector<C> &c,
         std::vector<D> &d, std::vector<std::string> &s, size_t n,
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
        A a;
        B b;
        C c;
        D d;
        string s;

        if (is >> a >> b >> c >> d and getline(is, s, eol)) {
            istringstream ist(s);

            ist >> s;
            if (!ist)
                s.erase();

            ta.push_back(a);
            tb.push_back(b);
            tc.push_back(c);
            td.push_back(d);
            ts.push_back(s);
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

#endif // RQ_READLINE_H
