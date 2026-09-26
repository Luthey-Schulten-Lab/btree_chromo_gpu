// Author(s): Ron Acda (using an iterative LLM-guided workflow, https://github.com/quarkron/iterative-hillclimber/tree/main)
// exact fast path for the text body of a custom dump (DumpCustom::convert_string, and write_lines when unbuffered).
//
// The hook appends the BD's final frame to <output>.lammpstrj with write_dump (id type x y z c_id_track c_type_track, 149 638
// atoms at t = 2996 s); LAMMPS formats it with one snprintf per field into the dump's string buffer: 0.20 s per hook. This
// file defines convert_string and write_lines (exported from the executable; DumpCustom calls them through member-function
// pointers, which bind to these definitions). When every column uses one of LAMMPS' default formats, "%d " (INT) or "%g "
// (DOUBLE, i.e. "%.6g "; the last column without the space), the fields are formatted with std::to_chars (decimal; general with precision 6, byte-identical to
// "%.6g", the to-chars change) into the same buffer, grown by the same per-line rule, with the same returned length and terminator; any
// other column layout goes to LAMMPS' own functions (RTLD_NEXT). WCM_FAST_DUMP_OFF=1: always LAMMPS. WCM_FAST_DUMP_VERIFY=1: each call
// is also rendered by LAMMPS' function and compared byte for byte (the fast text is what gets written).
#define protected public   // DumpCustom's column tables and Dump::fp are protected; the class layout is unchanged
#include "dump_custom.h"
#include "memory.h"
#undef protected
#include <climits>
#include <charconv>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <dlfcn.h>
#include <string>

namespace {

using WriteLines = void (*)(LAMMPS_NS::DumpCustom *, int, double *);

WriteLines real_write_lines()
{
  static WriteLines f = nullptr;
  if (!f) {
    void *p = dlsym(RTLD_NEXT, "_ZN9LAMMPS_NS10DumpCustom11write_linesEiPd");
    if (!p) { fprintf(stderr, "LAMMPS DumpCustom::write_lines not found (%s)\n", dlerror()); abort(); }
    std::memcpy(&f, &p, sizeof(p));
  }
  return f;
}

inline void put_int(std::string &out, int v)
{
  char b[16]; const std::to_chars_result r = std::to_chars(b, b + sizeof(b), v); out.append(b, (size_t)(r.ptr - b));
}

inline void put_g(std::string &out, double v)
{
  char b[40]; const std::to_chars_result r = std::to_chars(b, b + sizeof(b), v, std::chars_format::general, 6);
  out.append(b, (size_t)(r.ptr - b));
}

// the column layout the fast paths reproduce: every field "%d" (INT) or "%g" (DOUBLE), each followed by one space or by
// nothing (LAMMPS leaves the last column's format without the space)
bool plain_columns(const LAMMPS_NS::DumpCustom *d, bool *is_int, bool *space)
{
  const int nf = d->nfield;
  if (nf <= 0 || nf > 64) return false;
  for (int j = 0; j < nf; j++) {
    const char *f = d->vformat[j];
    if (!f) return false;
    space[j] = std::strcmp(f, "%d ") == 0 || std::strcmp(f, "%g ") == 0;
    const bool bare = std::strcmp(f, "%d") == 0 || std::strcmp(f, "%g") == 0;
    if (!space[j] && !bare) return false;
    if (d->vtype[j] == LAMMPS_NS::Dump::INT && f[1] == 'd') is_int[j] = true;
    else if (d->vtype[j] == LAMMPS_NS::Dump::DOUBLE && f[1] == 'g') is_int[j] = false;
    else return false;
  }
  return true;
}

inline int put_int_at(char *p, int v) { return (int)(std::to_chars(p, p + 16, v).ptr - p); }
inline int put_g_at(char *p, double v) { return (int)(std::to_chars(p, p + 40, v, std::chars_format::general, 6).ptr - p); }

}  // namespace

extern "C" int wcm_convert_string(LAMMPS_NS::DumpCustom *d, int n, double *mybuf) __asm__("_ZN9LAMMPS_NS10DumpCustom14convert_stringEiPd");
extern "C" int wcm_convert_string(LAMMPS_NS::DumpCustom *d, int n, double *mybuf)
{
  using F = int (*)(LAMMPS_NS::DumpCustom *, int, double *);
  static F real = nullptr;
  if (!real) {
    void *p = dlsym(RTLD_NEXT, "_ZN9LAMMPS_NS10DumpCustom14convert_stringEiPd");
    if (!p) { fprintf(stderr, "LAMMPS DumpCustom::convert_string not found (%s)\n", dlerror()); abort(); }
    std::memcpy(&real, &p, sizeof(p));
  }
  static const bool off = std::getenv("WCM_FAST_DUMP_OFF") != nullptr;
  static const bool verify = std::getenv("WCM_FAST_DUMP_VERIFY") != nullptr;
  bool is_int[64], space[64];
  if (off || !plain_columns(d, is_int, space)) return real(d, n, mybuf);
  std::string ref;
  if (verify) { const int r = real(d, n, mybuf); if (r < 0) return r; ref.assign(d->sbuf, (size_t)r); }
  constexpr int ONEFIELD = 32, DELTA = 1048576;   // dump_custom.cpp's constants: each field fits in ONEFIELD characters
  const int nf = d->nfield;
  int offset = 0, m = 0;
  for (int i = 0; i < n; i++) {
    if (offset + nf * ONEFIELD > d->maxsbuf) {
      if ((LAMMPS_NS::bigint) d->maxsbuf + DELTA > INT_MAX) return -1;
      d->maxsbuf += DELTA;
      d->memory->grow(d->sbuf, d->maxsbuf, "dump:sbuf");
    }
    char *p = d->sbuf;
    for (int j = 0; j < nf; j++) {
      offset += is_int[j] ? put_int_at(p + offset, static_cast<int>(mybuf[m])) : put_g_at(p + offset, mybuf[m]);
      if (space[j]) p[offset++] = ' ';
      m++;
    }
    p[offset++] = '\n';
  }
  d->sbuf[offset] = '\0';   // snprintf's terminator after the last "\n"
  if (verify)
    fprintf(stderr, "WCM_FAST_DUMP_VERIFY: convert_string %d lines %s (%d / %zu bytes)\n", n,
            (ref.size() == (size_t)offset && std::memcmp(ref.data(), d->sbuf, (size_t)offset) == 0) ? "IDENTICAL" : "DIFFERENT", offset, ref.size());
  return offset;
}

extern "C" void wcm_write_lines(LAMMPS_NS::DumpCustom *d, int n, double *mybuf) __asm__("_ZN9LAMMPS_NS10DumpCustom11write_linesEiPd");
extern "C" void wcm_write_lines(LAMMPS_NS::DumpCustom *d, int n, double *mybuf)
{
  static const bool off = std::getenv("WCM_FAST_DUMP_OFF") != nullptr;
  static const bool verify = std::getenv("WCM_FAST_DUMP_VERIFY") != nullptr;
  const int nf = d->nfield;
  bool is_int[64], space[64];
  const bool plain = !off && plain_columns(d, is_int, space);
  if (!plain) { real_write_lines()(d, n, mybuf); return; }

  std::string out;
  out.reserve((size_t)n * (size_t)nf * 12 + 16);
  int m = 0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < nf; j++) {
      if (is_int[j]) put_int(out, static_cast<int>(mybuf[m])); else put_g(out, mybuf[m]);
      if (space[j]) out.push_back(' ');
      m++;
    }
    out.push_back('\n');
  }
  if (verify) {
    FILE *keep = d->fp, *tmp = std::tmpfile();
    d->fp = tmp; real_write_lines()(d, n, mybuf); d->fp = keep;
    std::fflush(tmp); const long sz = std::ftell(tmp); std::rewind(tmp);
    std::string ref((size_t)(sz > 0 ? sz : 0), '\0');
    const size_t got = sz > 0 ? std::fread(&ref[0], 1, (size_t)sz, tmp) : 0;
    std::fclose(tmp);
    fprintf(stderr, "WCM_FAST_DUMP_VERIFY: write_lines %d lines %s (%zu / %zu bytes)\n", n, (got == ref.size() && ref == out) ? "IDENTICAL" : "DIFFERENT",
            out.size(), ref.size());
  }
  std::fwrite(out.data(), 1, out.size(), d->fp);
}
