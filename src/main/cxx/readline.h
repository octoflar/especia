// Input procedures
// Copyright (c) 2016, Ralf Quast
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
// TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
#ifndef RQ_READLINE_H
#define RQ_READLINE_H

#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace RQ {
    template<class A>
    std::istream& read(std::istream& is, std::vector<A>& a, size_t n,
        bool append = false);

    template<class A, class B>
    std::istream& read(std::istream& is, std::vector<A>& a, std::vector<B>& b, size_t n,
        bool append = false);

    template<class A, class B, class C>
    std::istream& read(std::istream& is, std::vector<A>& a, std::vector<B>& b,
        std::vector<C>& c, size_t n, bool append = false);

    template<class A, class B, class C, class D>
    std::istream& read(std::istream& is, std::vector<A>& a, std::vector<B>& b,
        std::vector<C>& c, std::vector<D>& d, size_t n, bool append = false);

    template<class A, class B, class C, class D, class E>
    std::istream& read(std::istream& is, std::vector<A>& a, std::vector<B>& b,
        std::vector<C>& c, std::vector<D>& d, std::vector<E>& e, size_t n,
        bool append = false);

    template<class A, class B, class C, class D>
    std::istream& read(std::istream& is, std::vector<A>& a, std::vector<B>& b,
        std::vector<C>& c, std::vector<D>& d, std::vector<std::string>& s, size_t n,
        char eol, bool append = false);

    // Read a line while skipping empty lines and stripping comments
    std::istream& readline(std::istream& is, std::string& s,
        char comment_mark = 0, char eol = '\n');
}

template<class A>
std::istream&
RQ::read(std::istream& is, std::vector<A>& a, size_t n, bool append)
{
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
std::istream&
RQ::read(std::istream& is, std::vector<A>& a, std::vector<B>& b, size_t n,
    bool append)
{
    using namespace std;

    vector<A> ta;
    vector<B> tb;

    ta.reserve(n);
    tb.reserve(n);

    for (size_t i = 0; i < n; ++i) {
        A a; B b;

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
std::istream&
RQ::read(std::istream& is, std::vector<A>& a, std::vector<B>& b, std::vector<C>& c,
    size_t n, bool append)
{
    using namespace std;

    vector<A> ta;
    vector<B> tb;
    vector<C> tc;

    ta.reserve(n);
    tb.reserve(n);
    tc.reserve(n);
    
    for (size_t i = 0; i < n; ++i) {
        A a; B b; C c;

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
std::istream&
RQ::read(std::istream& is, std::vector<A>& a, std::vector<B>& b, std::vector<C>& c,
    std::vector<D>& d, size_t n, bool append)
{
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
        A a; B b; C c; D d;

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
std::istream&
RQ::read(std::istream& is, std::vector<A>& a, std::vector<B>& b, std::vector<C>& c,
    std::vector<D>& d, std::vector<E>& e, size_t n, bool append)
{
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
        A a; B b; C c; D d; E e;

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
std::istream&
RQ::read(std::istream& is, std::vector<A>& a, std::vector<B>& b, std::vector<C>& c,
    std::vector<D>& d, std::vector<std::string>& s, size_t n,
    char eol, bool append)
{
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
        A a; B b; C c; D d; string s;

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
