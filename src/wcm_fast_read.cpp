// Author(s): Ron Acda (using an iterative LLM-guided workflow, https://github.com/quarkron/iterative-hillclimber/tree/main)
// exact fast path for the Bonds and Angles sections of read_data.
//
// read_data parses every Bonds and Angles line twice (a first pass that only counts entries per atom, then the store pass), and
// LAMMPS' Atom::data_bonds / Atom::data_angles tokenise each line into a vector of std::string and convert through
// utils::tnumeric / inumeric: 0.50 s of the hook's 0.63 s read at t = 2996 s (108 676 bonds, 108 676 angles). This file defines
// both functions (exported from the executable, so read_data binds to them). A chunk is taken by the fast path only when every
// line in it is plain: n newline-terminated lines of exactly 4 (5) whitespace-separated words, no comment, numeric type and atom
// IDs, no type labels, and every ID/type check LAMMPS makes passes. Those lines are then applied in order with the same
// arithmetic, the same map() lookups and the same stores (and data_bonds_post) as LAMMPS, so the atom arrays end up identical.
// Any other chunk goes unchanged to LAMMPS' own function (reached with RTLD_NEXT), which also raises LAMMPS' own errors.
// WCM_FAST_READ_BONDS_OFF=1: always LAMMPS' functions.
#include "atom.h"
#include "atom_vec.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <dlfcn.h>
#include <vector>

using LAMMPS_NS::tagint;

// force->newton_bond of the instance whose read_data is running (Atom::lmp is protected); set by LAMMPS_simulator::read_data
// around its read_data, -1 otherwise (then every chunk goes to LAMMPS)
int wcm_read_newton_bond = -1;

namespace {

bool fast_off() { static const bool off = std::getenv("WCM_FAST_READ_BONDS_OFF") != nullptr; return off; }

template <class F> F real_fn(const char *sym)
{
  void *p = dlsym(RTLD_NEXT, sym);
  if (!p) { fprintf(stderr, "LAMMPS %s not found (%s)\n", sym, dlerror()); abort(); }
  F f; std::memcpy(&f, &p, sizeof(p)); return f;
}

// the separators of LAMMPS' Tokenizer (" \t\r\n\f"; the newline ends the line here)
inline bool is_space(char c) { return c == ' ' || c == '\t' || c == '\r' || c == '\f'; }

// ^[+-]?\d+$ as utils::is_integer, converted as strtoll (tnumeric) / strtol (inumeric) do for such a token
bool parse_int(const char *b, const char *e, long long &v)
{
  const char *p = b;
  if (p < e && (*p == '+' || *p == '-')) p++;
  if (p == e) return false;
  for (const char *q = p; q < e; q++) if (*q < '0' || *q > '9') return false;
  if (e - p > 18) return false;          // let LAMMPS handle anything that could overflow
  long long x = 0;
  for (const char *q = p; q < e; q++) x = x * 10 + (*q - '0');
  v = (*b == '-') ? -x : x;
  return true;
}

// Splits n newline-terminated lines of buf into words; fills rows of W integers (words 1..W-1 converted, word 0 ignored as
// LAMMPS ignores it). Returns false (nothing touched) if any line is not plain: wrong word count, a comment, a non-integer.
// Blank lines are skipped as LAMMPS skips lines with zero words.
template <int W>
bool split_chunk(int n, const char *buf, std::vector<long long> &rows)
{
  rows.clear(); rows.reserve((size_t)n * (W - 1));
  const char *p = buf;
  for (int i = 0; i < n; i++) {
    const char *nl = std::strchr(p, '\n');
    if (!nl) return false;
    const char *wb[W + 1], *we[W + 1]; int nw = 0;
    const char *q = p;
    while (q < nl) {
      while (q < nl && is_space(*q)) q++;
      if (q >= nl) break;
      const char *s = q;
      while (q < nl && !is_space(*q)) q++;
      if (*s == '#') break;               // a comment ends the words (the fast path then requires the count to be W)
      if (nw == W) return false;
      wb[nw] = s; we[nw] = q; nw++;
    }
    for (const char *c = p; c < nl; c++) if (*c == '#') return false;   // any comment at all: LAMMPS decides
    if (nw != 0) {
      if (nw != W) return false;
      for (int w = 1; w < W; w++) { long long v; if (!parse_int(wb[w], we[w], v)) return false; rows.push_back(v); }
    }
    p = nl + 1;
  }
  return true;
}

}  // namespace

extern "C" void wcm_data_bonds(LAMMPS_NS::Atom *a, int n, char *buf, int *count, tagint id_offset, int type_offset,
                                int labelflag, int *ilabel) __asm__("_ZN9LAMMPS_NS4Atom10data_bondsEiPcPiiiiS2_");
extern "C" void wcm_data_bonds(LAMMPS_NS::Atom *a, int n, char *buf, int *count, tagint id_offset, int type_offset,
                                int labelflag, int *ilabel)
{
  using F = void (*)(LAMMPS_NS::Atom *, int, char *, int *, tagint, int, int, int *);
  static F real = real_fn<F>("_ZN9LAMMPS_NS4Atom10data_bondsEiPcPiiiiS2_");
  static std::vector<long long> rows;
  if (fast_off() || wcm_read_newton_bond < 0 || labelflag || !split_chunk<4>(n, buf, rows)) { real(a, n, buf, count, id_offset, type_offset, labelflag, ilabel); return; }
  const size_t nl = rows.size() / 3;
  for (size_t k = 0; k < nl; k++) {       // every check LAMMPS makes, before anything is stored
    long long t = rows[3 * k] + type_offset, a1 = rows[3 * k + 1] + id_offset, a2 = rows[3 * k + 2] + id_offset;
    if (t < 1 || t > a->nbondtypes || a1 <= 0 || a1 > a->map_tag_max || a2 <= 0 || a2 > a->map_tag_max || a1 == a2)
      { real(a, n, buf, count, id_offset, type_offset, labelflag, ilabel); return; }
  }
  const int newton_bond = wcm_read_newton_bond;
  for (size_t k = 0; k < nl; k++) {
    const int itype = (int)(rows[3 * k] + type_offset);
    const tagint atom1 = (tagint)(rows[3 * k + 1] + id_offset), atom2 = (tagint)(rows[3 * k + 2] + id_offset);
    int m;
    if ((m = a->map(atom1)) >= 0) {
      if (count) count[m]++;
      else {
        a->bond_type[m][a->num_bond[m]] = itype; a->bond_atom[m][a->num_bond[m]] = atom2; a->num_bond[m]++;
        a->avec->data_bonds_post(m, a->num_bond[m], atom1, atom2, id_offset);
      }
    }
    if (newton_bond == 0 && (m = a->map(atom2)) >= 0) {
      if (count) count[m]++;
      else {
        a->bond_type[m][a->num_bond[m]] = itype; a->bond_atom[m][a->num_bond[m]] = atom1; a->num_bond[m]++;
        a->avec->data_bonds_post(m, a->num_bond[m], atom1, atom2, id_offset);
      }
    }
  }
}

extern "C" void wcm_data_angles(LAMMPS_NS::Atom *a, int n, char *buf, int *count, tagint id_offset, int type_offset,
                                 int labelflag, int *ilabel) __asm__("_ZN9LAMMPS_NS4Atom11data_anglesEiPcPiiiiS2_");
extern "C" void wcm_data_angles(LAMMPS_NS::Atom *a, int n, char *buf, int *count, tagint id_offset, int type_offset,
                                 int labelflag, int *ilabel)
{
  using F = void (*)(LAMMPS_NS::Atom *, int, char *, int *, tagint, int, int, int *);
  static F real = real_fn<F>("_ZN9LAMMPS_NS4Atom11data_anglesEiPcPiiiiS2_");
  static std::vector<long long> rows;
  if (fast_off() || wcm_read_newton_bond < 0 || labelflag || !split_chunk<5>(n, buf, rows)) { real(a, n, buf, count, id_offset, type_offset, labelflag, ilabel); return; }
  const size_t nl = rows.size() / 4;
  for (size_t k = 0; k < nl; k++) {
    long long t = rows[4 * k] + type_offset, a1 = rows[4 * k + 1] + id_offset, a2 = rows[4 * k + 2] + id_offset, a3 = rows[4 * k + 3] + id_offset;
    if (t < 1 || t > a->nangletypes || a1 <= 0 || a1 > a->map_tag_max || a2 <= 0 || a2 > a->map_tag_max ||
        a3 <= 0 || a3 > a->map_tag_max || a1 == a2 || a1 == a3 || a2 == a3)
      { real(a, n, buf, count, id_offset, type_offset, labelflag, ilabel); return; }
  }
  const int newton_bond = wcm_read_newton_bond;
  auto store = [a](int m, int itype, tagint a1, tagint a2, tagint a3) {
    a->angle_type[m][a->num_angle[m]] = itype; a->angle_atom1[m][a->num_angle[m]] = a1;
    a->angle_atom2[m][a->num_angle[m]] = a2; a->angle_atom3[m][a->num_angle[m]] = a3; a->num_angle[m]++;
  };
  for (size_t k = 0; k < nl; k++) {
    const int itype = (int)(rows[4 * k] + type_offset);
    const tagint a1 = (tagint)(rows[4 * k + 1] + id_offset), a2 = (tagint)(rows[4 * k + 2] + id_offset), a3 = (tagint)(rows[4 * k + 3] + id_offset);
    int m;
    if ((m = a->map(a2)) >= 0) { if (count) count[m]++; else store(m, itype, a1, a2, a3); }
    if (newton_bond == 0) {
      if ((m = a->map(a1)) >= 0) { if (count) count[m]++; else store(m, itype, a1, a2, a3); }
      if ((m = a->map(a3)) >= 0) { if (count) count[m]++; else store(m, itype, a1, a2, a3); }
    }
  }
}

// ============================================================================================================================
// exact fast path for the Atoms section of read_data (Atom::data_atoms, 0.21 s of the 0.33 s read at t = 2996 s,
// 149 638 atoms). LAMMPS tokenises every line into a vector of std::string (twice for the first line), matches comments with
// utils::strmatch and converts each value through utils::numeric / inumeric (std::string copies + regex-style checks +
// std::stod / std::stoi). A chunk is taken by the fast path only when all of these hold: a 3d orthogonal box that is
// non-periodic in x, y and z (Domain::remap is then the identity and no sub-box epsilon applies), no general triclinic, no
// shift, no type labels, atom IDs enabled, an atom style whose Atoms line is "id molecule type x y z" without image flags
// (size_data_atom 6, x in column 4, molecule flag, 4-byte tagint), and every line of the chunk plain: newline-terminated,
// exactly 6 words, no comment, integer id / molecule within int, an all-digit type within 1..ntypes, coordinates that match
// utils::is_double and convert without range error, a positive id. Those lines are then applied in order with exactly what
// Atom::data_atoms + AtomVec::data_atom do for them (owned-atom test against the sub-box, grow at nmax, x, mask 1, image of
// flags 0, v 0, tag, molecule, the style's data_atom_post, nlocal++, id/molecule offsets, type). A coordinate is the
// correctly rounded value of its decimal token, as std::stod (strtod) gives it: by Clinger's exact shortcut when the token has
// at most 15 significant digits and a decimal exponent within +-22 (all of btree's %.6g text), by strtod otherwise. Any other
// chunk goes to LAMMPS' own function (RTLD_NEXT). WCM_FAST_READ_ATOMS_VERIFY=1 re-reads with LAMMPS' parser and compares the LAMMPS
// state byte for byte (btree_driver::simulator_read_data). WCM_FAST_READ_ATOMS_OFF=1: always LAMMPS' function.
// ============================================================================================================================
#include "lammps.h"
#include "domain.h"
#include "comm.h"
#include <cerrno>
#include <cstring>

LAMMPS_NS::LAMMPS *wcm_read_lmp = nullptr;   // the instance whose read_data is running (set by LAMMPS_simulator::read_data)
bool wcm_fast_read_atoms_force_off = false;               // WCM_FAST_READ_ATOMS_VERIFY's reference read (LAMMPS' own parser)

namespace {

bool fast_fast_read_atoms_off() { static const bool off = std::getenv("WCM_FAST_READ_ATOMS_OFF") != nullptr; return off; }

// ^[+-]?\d+$ within int (utils::is_integer + std::stoi)
bool parse_int32(const char *b, const char *e, int &v)
{
  long long x;
  if (!parse_int(b, e, x) || x < -2147483648LL || x > 2147483647LL) return false;
  v = (int)x; return true;
}

// utils::is_double: ^[+-]?(\d+\.?\d*|\d*\.?\d+)([eE][+-]?\d+)?$, then std::stod (strtod over the whole token, ERANGE = error)
bool parse_double(const char *b, const char *e, double &v)
{
  const char *p = b;
  if (p < e && (*p == '+' || *p == '-')) p++;
  int nd = 0; bool dot = false;
  while (p < e && ((*p >= '0' && *p <= '9') || (*p == '.' && !dot))) { if (*p == '.') dot = true; else nd++; p++; }
  if (nd == 0) return false;
  if (p < e) {
    if (*p != 'e' && *p != 'E') return false;
    p++;
    if (p < e && (*p == '+' || *p == '-')) p++;
    const char *d = p;
    while (p < e && *p >= '0' && *p <= '9') p++;
    if (p == d || p != e) return false;
  }
  // Exact shortcut (Clinger): at most 15 significant digits make an integer M < 2^53 held exactly; with a decimal exponent
  // |k| <= 22, 10^|k| is exact too, so M * 10^k or M / 10^-k is one correctly rounded IEEE operation = the correctly rounded
  // value of the decimal, which is what strtod returns. Anything else goes to strtod.
  {
    static const double p10[23] = {1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12, 1e13, 1e14, 1e15,
                                   1e16, 1e17, 1e18, 1e19, 1e20, 1e21, 1e22};
    const char *q = b; bool neg = false;
    if (*q == '+' || *q == '-') { neg = (*q == '-'); q++; }
    unsigned long long m = 0; int sig = 0, frac = 0; bool in_frac = false, lead = true;
    for (; q < e && *q != 'e' && *q != 'E'; q++) {
      if (*q == '.') { in_frac = true; continue; }
      const int dgt = *q - '0';
      if (in_frac) frac++;
      if (lead && dgt == 0) continue;     // leading zeros are not significant
      lead = false;
      if (++sig > 15) break;
      m = m * 10 + (unsigned)dgt;
    }
    if (sig <= 15) {
      long long ex = 0;
      if (q < e) {                        // exponent (syntax checked above)
        q++; bool eneg = false;
        if (*q == '+' || *q == '-') { eneg = (*q == '-'); q++; }
        for (; q < e && ex < 100000; q++) ex = ex * 10 + (*q - '0');
        if (eneg) ex = -ex;
      }
      const long long k = ex - frac;
      if (k >= -22 && k <= 22) {
        const double md = (double)m;
        const double r = (k >= 0) ? md * p10[k] : md / p10[-k];
        v = neg ? -r : r;
        return true;
      }
    }
  }
  char tmp[64];
  if (e - b >= (long)sizeof(tmp)) return false;
  std::memcpy(tmp, b, e - b); tmp[e - b] = '\0';
  char *end = nullptr;
  errno = 0;
  v = std::strtod(tmp, &end);
  return errno != ERANGE && end == tmp + (e - b);
}

struct AtomRow { int id, mol, type; double x[3]; };

// n newline-terminated lines of plain "id molecule type x y z" (blank lines skipped, as LAMMPS skips lines with zero words;
// the first line must be a data line, as LAMMPS takes the word count from it)
bool split_atoms(int n, const char *buf, std::vector<AtomRow> &rows)
{
  rows.clear(); rows.reserve(n);
  const char *p = buf;
  for (int i = 0; i < n; i++) {
    const char *nl = std::strchr(p, '\n');
    if (!nl) return false;
    const char *wb[6], *we[6]; int nw = 0;
    const char *q = p;
    while (q < nl) {
      while (q < nl && is_space(*q)) q++;
      if (q >= nl) break;
      const char *s = q;
      while (q < nl && !is_space(*q)) q++;
      if (nw == 6) return false;
      wb[nw] = s; we[nw] = q; nw++;
    }
    if (nw == 0) { if (i == 0) return false; p = nl + 1; continue; }
    if (nw != 6) return false;
    AtomRow r;
    if (!parse_int32(wb[0], we[0], r.id) || !parse_int32(wb[1], we[1], r.mol)) return false;
    for (const char *c = wb[2]; c < we[2]; c++) if (*c < '0' || *c > '9') return false;
    if (!parse_int32(wb[2], we[2], r.type)) return false;
    for (int d = 0; d < 3; d++) if (!parse_double(wb[3 + d], we[3 + d], r.x[d])) return false;
    rows.push_back(r);
    p = nl + 1;
  }
  return true;
}

}  // namespace

extern "C" void wcm_data_atoms(LAMMPS_NS::Atom *a, int n, char *buf, tagint id_offset, tagint mol_offset, int type_offset,
                                int shiftflag, double *shift, int labelflag, int *ilabel, int triclinic_general)
                                __asm__("_ZN9LAMMPS_NS4Atom10data_atomsEiPciiiiPdiPii");
extern "C" void wcm_data_atoms(LAMMPS_NS::Atom *a, int n, char *buf, tagint id_offset, tagint mol_offset, int type_offset,
                                int shiftflag, double *shift, int labelflag, int *ilabel, int triclinic_general)
{
  using F = void (*)(LAMMPS_NS::Atom *, int, char *, tagint, tagint, int, int, double *, int, int *, int);
  static F real = real_fn<F>("_ZN9LAMMPS_NS4Atom10data_atomsEiPciiiiPdiPii");
  static std::vector<AtomRow> rows;
  LAMMPS_NS::LAMMPS *lmp = wcm_read_lmp;
  auto fallback = [&]() { real(a, n, buf, id_offset, mol_offset, type_offset, shiftflag, shift, labelflag, ilabel, triclinic_general); };
  if (fast_fast_read_atoms_off() || wcm_fast_read_atoms_force_off || !lmp || lmp->atom != a || labelflag || shiftflag || triclinic_general) { fallback(); return; }
  LAMMPS_NS::Domain *dom = lmp->domain;
  if (sizeof(tagint) != 4 || dom->dimension != 3 || dom->triclinic || dom->xperiodic || dom->yperiodic || dom->zperiodic ||
      !a->tag_enable || !a->molecule_flag || !a->avec || a->avec->size_data_atom != 6 || a->avec->xcol_data != 4 ||
      a->labelmapflag)
    { fallback(); return; }
  if (!split_atoms(n, buf, rows)) { fallback(); return; }
  for (const AtomRow &r : rows) {         // every check LAMMPS makes, before anything is stored
    const long long t = (long long)r.type + type_offset;
    if (r.id <= 0 || t < 1 || t > a->ntypes) { fallback(); return; }
  }
  const double *sublo = dom->sublo, *subhi = dom->subhi;
  const LAMMPS_NS::imageint imagedata = ((LAMMPS_NS::imageint)IMGMAX & IMGMASK) | (((LAMMPS_NS::imageint)IMGMAX & IMGMASK) << IMGBITS) |
                                        (((LAMMPS_NS::imageint)IMGMAX & IMGMASK) << IMG2BITS);
  for (const AtomRow &r : rows) {
    const double *c = r.x;
    if (!(c[0] >= sublo[0] && c[0] < subhi[0] && c[1] >= sublo[1] && c[1] < subhi[1] && c[2] >= sublo[2] && c[2] < subhi[2])) continue;
    const int nlocal = a->nlocal;
    if (nlocal == a->nmax) a->avec->grow(0);
    a->x[nlocal][0] = c[0]; a->x[nlocal][1] = c[1]; a->x[nlocal][2] = c[2];
    a->mask[nlocal] = 1;
    a->image[nlocal] = imagedata;
    a->v[nlocal][0] = 0.0; a->v[nlocal][1] = 0.0; a->v[nlocal][2] = 0.0;
    a->tag[nlocal] = r.id;
    a->molecule[nlocal] = r.mol;
    a->avec->data_atom_post(nlocal);
    a->nlocal++;
    if (id_offset) a->tag[nlocal] += id_offset;
    if (mol_offset) a->molecule[nlocal] += mol_offset;
    a->type[nlocal] = r.type + type_offset;
  }
}
