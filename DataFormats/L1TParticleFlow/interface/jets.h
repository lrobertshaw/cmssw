#ifndef DataFormats_L1TParticleFlow_jets_h
#define DataFormats_L1TParticleFlow_jets_h

#include "DataFormats/L1TParticleFlow/interface/datatypes.h"
#include "DataFormats/L1TParticleFlow/interface/gt_datatypes.h"
#include "DataFormats/L1TParticleFlow/interface/bit_encoding.h"
#include <array>
#include <cstdint>

namespace l1ct {

  struct Jet {
    pt_t hwPt;
    glbeta_t hwEta;
    glbphi_t hwPhi;
    z0_t hwZ0;
    // b_tag_score_t hwBtagScore;
    mass_t hwMass;

    inline bool operator==(const Jet &other) const {
      return hwPt == other.hwPt && hwEta == other.hwEta && hwPhi == other.hwPhi;
    }

    inline bool operator>(const Jet &other) const { return hwPt > other.hwPt; }
    inline bool operator<(const Jet &other) const { return hwPt < other.hwPt; }

    inline void clear() {
      hwPt = 0;
      hwEta = 0;
      hwPhi = 0;
      hwZ0 = 0;
      // hwBtagScore = 0;
      hwMass = 0;
    }

    int intPt() const { return Scales::intPt(hwPt); }
    int intEta() const { return hwEta.to_int(); }
    int intPhi() const { return hwPhi.to_int(); }
    float floatPt() const { return Scales::floatPt(hwPt); }
    float floatEta() const { return Scales::floatEta(hwEta); }
    float floatPhi() const { return Scales::floatPhi(hwPhi); }
    float floatZ0() const { return Scales::floatZ0(hwZ0); }
    // float floatBtagScore() const { return Scales::floatBtagScore(hwBtagScore); }
    float floatMass() const { return Scales::floatMass(hwMass); }

    static const int BITWIDTH = pt_t::width + glbeta_t::width + glbphi_t::width + z0_t::width + mass_t::width; //b_tag_score_t::width;    //l1ct types give 59, but vhdl expects 57
    inline ap_uint<BITWIDTH> pack_ap() const {
      ap_uint<BITWIDTH> ret;
      unsigned int start = 0;
      pack_into_bits(ret, start, hwPt);
      pack_into_bits(ret, start, hwEta);
      pack_into_bits(ret, start, hwPhi);
      pack_into_bits(ret, start, hwZ0);
      // pack_into_bits(ret, start, hwBtagScore);
      pack_into_bits(ret, start, hwMass);
      return ret;
    }

    inline std::array<uint64_t, 2> pack() const {
      std::array<uint64_t, 2> packed = {{0, 0}};
      ap_uint<BITWIDTH> bits = this->pack_ap();
      packed[0] = bits;
      //packed[1] = bits[slice]; // for when there are more than 64 bits in the word
      return packed;
    }

    inline static Jet unpack_ap(const ap_uint<BITWIDTH> &src) {
      Jet ret;
      ret.initFromBits(src);
      return ret;
    }

    inline void initFromBits(const ap_uint<BITWIDTH> &src) {
      unsigned int start = 0;
      unpack_from_bits(src, start, hwPt);
      unpack_from_bits(src, start, hwEta);
      unpack_from_bits(src, start, hwPhi);
      unpack_from_bits(src, start, hwZ0);
      // unpack_from_bits(src, start, hwBtagScore);
      unpack_from_bits(src, start, hwMass);
    }

    inline static Jet unpack(const std::array<uint64_t, 2> &src) {
      // just one set while the word has fewer than 64 bits
      ap_uint<BITWIDTH> bits = src[0];
      return unpack_ap(bits);
    }

    inline static Jet unpack(long long unsigned int &src) {
      // unpack from single 64b int
      ap_uint<BITWIDTH> bits = src;
      return unpack_ap(bits);
    }

    l1gt::Jet toGT() const {
      l1gt::Jet j;
      j.valid = hwPt != 0;
      j.v3.pt = CTtoGT_pt(hwPt);
      j.v3.phi = CTtoGT_phi(hwPhi);
      j.v3.eta = CTtoGT_eta(hwEta);
      j.z0(l1ct::z0_t::width - 1, 0) = hwZ0(l1ct::z0_t::width - 1, 0);
      // j.hwBtagScore = hwBtagScore;
      j.hwMass = CTtoGT_mass(hwMass);
      return j;
    }
  };

  inline void clear(Jet &c) { c.clear(); }

}  // namespace l1ct

#endif
