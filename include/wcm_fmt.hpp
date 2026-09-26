// Author(s): Ron Acda (using an iterative LLM-guided workflow, https://github.com/quarkron/iterative-hillclimber/tree/main)
#pragma once
// text formatting identical to std::ostream's defaults (int: decimal; double: precision 6, no floatfield = "%.6g",
// which is the format libstdc++ passes to vsnprintf), appended to a string.
#include <string>
#include <charconv>
#include <cstdio>
#include <cstdlib>

inline void wcm_put_int(std::string &out, long long v)
{
    char b[24]; int n = 0; bool neg = v < 0; unsigned long long u = neg ? 0ULL - (unsigned long long)v : (unsigned long long)v;
    do { b[n++] = (char)('0' + (u % 10)); u /= 10; } while (u);
    if (neg) out.push_back('-');
    while (n) out.push_back(b[--n]);
}

// std::to_chars(general, 6) is specified as printf's "%.6g" in the C locale; libstdc++ 11 gives the same bytes on
// 6.36 M values (the t=2996 s coordinates, 6 M random doubles over every exponent, rounding and special-value edges) at a
// third of snprintf's time. WCM_TO_CHARS_OFF=1: snprintf.
inline void wcm_put_double(std::string &out, double v)
{
    static const bool off_to_chars = std::getenv("WCM_TO_CHARS_OFF") != nullptr;
    char b[40];
    if (!off_to_chars) { const std::to_chars_result r = std::to_chars(b, b + sizeof(b), v, std::chars_format::general, 6); out.append(b, (size_t)(r.ptr - b)); return; }
    const int n = std::snprintf(b, sizeof(b), "%.6g", v);
    out.append(b, n > 0 ? (size_t)n : 0);
}
