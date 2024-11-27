#include "L1Trigger/Phase2L1ParticleFlow/interface/jetmet/L1SeedConePFJetEmulator.h"

L1SCJetEmu::L1SCJetEmu(bool debug, float coneSize, unsigned nJets)
    : debug_(debug),
      coneSize_(coneSize),
      nJets_(nJets),
      rCone2_(coneSize * coneSize / l1ct::Scales::ETAPHI_LSB / l1ct::Scales::ETAPHI_LSB) {
  init_invert_table<pt_t, inv_pt_t, N_table_inv_pt>(inv_pt_table_);
}

L1SCJetEmu::detaphi_t L1SCJetEmu::deltaPhi(L1SCJetEmu::Particle a, L1SCJetEmu::Particle b) {
  detaphi_t dphi = detaphi_t(a.hwPhi) - detaphi_t(b.hwPhi);
  // phi wrap
  detaphi_t dphi0 =
      dphi > detaphi_t(l1ct::Scales::INTPHI_PI) ? detaphi_t(l1ct::Scales::INTPHI_TWOPI - dphi) : detaphi_t(dphi);
  detaphi_t dphi1 =
      dphi < detaphi_t(-l1ct::Scales::INTPHI_PI) ? detaphi_t(l1ct::Scales::INTPHI_TWOPI + dphi) : detaphi_t(dphi);
  detaphi_t dphiw = dphi > detaphi_t(0) ? dphi0 : dphi1;
  return dphiw;
}

bool L1SCJetEmu::inCone(L1SCJetEmu::Particle seed, L1SCJetEmu::Particle part) const {
  // scale the particle eta, phi to hardware units
  detaphi_t deta = detaphi_t(seed.hwEta) - detaphi_t(part.hwEta);
  detaphi_t dphi = deltaPhi(seed, part);
  bool ret = deta * deta + dphi * dphi < rCone2_;
  //bool ret = r2 < cone2;
  if (debug_) {
    detaphi2_t r2 = detaphi2_t(deta) * detaphi2_t(deta) + detaphi2_t(dphi) * detaphi2_t(dphi);
    dbgCout() << "  part eta, seed eta: " << part.hwEta << ", " << seed.hwEta << std::endl;
    dbgCout() << "  part phi, seed phi: " << part.hwPhi << ", " << seed.hwPhi << std::endl;
    dbgCout() << "  pt, deta, dphi, r2, cone2, lt: " << part.hwPt << ", " << deta << ", " << dphi << ", "
              << deta * deta + dphi * dphi << ", " << rCone2_ << ", " << ret << std::endl;
  }
  return ret;
}

L1SCJetEmu::Jet L1SCJetEmu::makeJet_HW(const std::vector<Particle>& parts) const {
  // Seed Cone Jet algorithm with ap_fixed types and hardware emulation
  Particle seed = reduce(parts, op_max);    // find maximum pt particle

  // Event with saturation, order of terms doesn't matter since they're all positive
  auto sumpt = [](pt_t(a), const Particle& b) { return a + b.hwPt; };    // essentially a python lambda fn

  // Sum the pt
  pt_t pt = std::accumulate(parts.begin(), parts.end(), pt_t(0), sumpt);
  inv_pt_t inv_pt = invert_with_shift<pt_t, inv_pt_t, N_table_inv_pt>(pt, inv_pt_table_, false);

  // pt weighted d eta
  std::vector<pt_etaphi_t> pt_deta;
  pt_deta.resize(parts.size());
  std::transform( parts.begin(), parts.end(), pt_deta.begin(), [&seed](const Particle& part) {
    // In the firmware we calculate the per-particle pt-weighted deta
    return pt_etaphi_t(part.hwPt * detaphi_t(part.hwEta - seed.hwEta));
  } );
  // Accumulate the pt-weighted etas. Init to 0, include seed in accumulation
  pt_etaphi_t sum_pt_eta = std::accumulate(pt_deta.begin(), pt_deta.end(), pt_etaphi_t(0));
  etaphi_t eta = seed.hwEta + etaphi_t(sum_pt_eta * inv_pt);

  // pt weighted d phi
  std::vector<pt_etaphi_t> pt_dphi;
  pt_dphi.resize(parts.size());
  std::transform(parts.begin(), parts.end(), pt_dphi.begin(), [&seed](const Particle& part) {
    // In the firmware we calculate the per-particle pt-weighted dphi
    return pt_etaphi_t(part.hwPt * deltaPhi(part, seed));
  });
  // Accumulate the pt-weighted phis. Init to 0, include seed in accumulation
  pt_etaphi_t sum_pt_phi = std::accumulate(pt_dphi.begin(), pt_dphi.end(), pt_etaphi_t(0));
  etaphi_t phi = seed.hwPhi + etaphi_t(sum_pt_phi * inv_pt);    // shift the seed by pt weighted sum_pt_phi

  std::vector<Particle> truncated;
  std::copy(parts.begin(), parts.end(), back_inserter(truncated));
  std::sort(truncated.begin(), truncated.end(), [](const Particle& a, const Particle& b) {return a.hwPt > b.hwPt;});    // sort by pt by jet mass fn as in firmware
  truncated.resize(NCONSTITS);    // truncate to 4 particles
  mass_t mass = L1SCJetEmu::jetMass_HW( truncated, seed );

  Jet jet;
  jet.hwPt = pt;
  jet.hwEta = eta;
  jet.hwPhi = phi;
  jet.hwMass = mass;
  jet.constituents = parts;

  if (debug_) {
    std::for_each(pt_dphi.begin(), pt_dphi.end(), [](pt_etaphi_t& x) { dbgCout() << "pt_dphi: " << x << std::endl; });
    std::for_each(pt_deta.begin(), pt_deta.end(), [](pt_etaphi_t& x) { dbgCout() << "pt_deta: " << x << std::endl; });
    dbgCout() << " sum_pt_eta: " << sum_pt_eta << ", 1/pt: " << inv_pt
              << ", sum_pt_eta * 1/pt: " << etaphi_t(sum_pt_eta * inv_pt) << std::endl;
    dbgCout() << " sum_pt_phi: " << sum_pt_phi << ", 1/pt: " << inv_pt
              << ", sum_pt_phi * 1/pt: " << etaphi_t(sum_pt_phi * inv_pt) << std::endl;
    dbgCout() << " uncorr eta: " << seed.hwEta << ", phi: " << seed.hwPhi << std::endl;
    dbgCout() << "   corr eta: " << eta << ", phi: " << phi << std::endl;
    dbgCout() << "         pt: " << pt << std::endl;
  }

  return jet;
}

L1SCJetEmu::mass_t L1SCJetEmu::jetMass_HW(const std::vector<Particle>& parts, const Particle& seed) const {    // need ampersand?

  // INSTANTIATE LUTS
  static constexpr int N = 185;
  static eventrig_t cosh_lut[N];
  static eventrig_t cos_lut[N];
  static oddtrig_t sin_lut[N];
  static oddtrig_t sinh_lut[N];
  for (unsigned hwEtaPhi = 0; hwEtaPhi < N; hwEtaPhi++) {
    float x = l1ct::Scales::floatEta((etaphi_t)hwEtaPhi);    // for each step in hardware units, convert
    cosh_lut[hwEtaPhi] = cosh(x); // Store cosh(hwEta) in hardware units
    cos_lut[hwEtaPhi] = cos(x);   // Store cos(hwEta) in hardware units
    sin_lut[hwEtaPhi] = sin(x);   // Store sin(hwEta) in hardware units
    sinh_lut[hwEtaPhi] = sinh(x); // Store sinh(hwEta) in hardware units
  }

  std::vector<ppt_t> energy;
  energy.resize(parts.size());
  std::transform(parts.begin(), parts.end(), energy.begin(), [&seed](const Particle& part) {
    detaphi_t etaIdx = part.hwEta - seed.hwEta;
    return ppt_t( part.hwPt * cosh_lut[std::abs(etaIdx)] );
  });
  ppt_t sum_energy = std::accumulate(energy.begin(), energy.end(), ppt_t(0));

  std::vector<ppt_t> px;
  px.resize(parts.size());
  std::transform(parts.begin(), parts.end(), px.begin(), [&seed](const Particle& part) {
    detaphi_t phiIdx = deltaPhi(part, seed);
    return ppt_t( part.hwPt * cos_lut[std::abs(phiIdx)] );
  });
  ppt_t sum_px = std::accumulate(px.begin(), px.end(), ppt_t(0));

  std::vector<npt_t> py;
  py.resize(parts.size());
  std::transform(parts.begin(), parts.end(), py.begin(), [&seed](const Particle& part) {
    detaphi_t phiIdx = deltaPhi(part, seed);
    int sign = (phiIdx >= 0) ? 1 : -1;
    return npt_t( part.hwPt * sin_lut[std::abs(phiIdx)] * sign );
  });
  npt_t sum_py = std::accumulate(py.begin(), py.end(), npt_t(0));

  std::vector<npt_t> pz;
  pz.resize(parts.size());
  std::transform(parts.begin(), parts.end(), pz.begin(), [&seed](const Particle& part) {
    detaphi_t etaIdx = part.hwEta - seed.hwEta;
    int sign = (etaIdx >= 0) ? 1 : -1;
    return npt_t( part.hwPt * sinh_lut[std::abs(etaIdx)] * sign );
  });
  npt_t sum_pz = std::accumulate(pz.begin(), pz.end(), npt_t(0));

  mass2_t mass2 = (sum_energy * sum_energy) - (sum_px * sum_px) - (sum_py * sum_py) - (sum_pz * sum_pz);
  mass_t mass = std::pow( mass2, 0.5 );
  
  return mass;
}


std::vector<L1SCJetEmu::Jet> L1SCJetEmu::emulateEvent(std::vector<Particle>& parts, std::vector<Particle>& seeds, bool useExternalSeeds, bool allowDoubleCounting) const {
  // The fixed point algorithm emulation
  std::vector<Particle> work;
  work.resize(parts.size());
  std::transform(parts.begin(), parts.end(), work.begin(), [](const Particle& part) { return part; });

  std::vector<Jet> jets;
  jets.reserve(nJets_);
  while ( ( !work.empty() && jets.size() < nJets_ )  ) {
    if ( useExternalSeeds && seeds.size() == 0 ) break; // Can this be combined with previous line?
    // Take the highest pt candidate as a seed
    // Use the firmware reduce function to find the same seed as the firmware
    // in case there are multiple seeds with the same pT
    // ... or use external seed if configured to do so
    Particle seed = (useExternalSeeds) ? seeds.at(0) : reduce(work, op_max);

    // Get the particles within a coneSize_ of the seed
    std::vector<Particle> particlesInCone;
    std::copy_if(work.begin(), work.end(), std::back_inserter(particlesInCone), [&](const Particle& part) {
      return inCone(seed, part);
    });
    if (debug_) {
      dbgCout() << "Seed: " << seed.hwPt << ", " << seed.hwEta << ", " << seed.hwPhi << std::endl;
      std::cout << "N particles : " << particlesInCone.size() << std::endl;
      std::for_each(particlesInCone.begin(), particlesInCone.end(), [&](Particle& part) {
        dbgCout() << "  Part: " << part.hwPt << ", " << part.hwEta << ", " << part.hwPhi << std::endl;
        inCone(seed, part);
      });
    }
    if ( particlesInCone.size() > 0 ) { // Possible hack - some seeds don't have any clustered particles.  Need to understand if real effect (could be) or a bug
     
    // if (debug_){ dbgCout() << "Internal seed used!" << std::endl;}   //debug
     jets.push_back(makeJet_HW(particlesInCone));
     //remove the clustered particles
     if(! allowDoubleCounting){       //if double couting not allowed
    //  if (debug_){ dbgCout() << "Removing candidates!" << std::endl;}   //debug
      work.erase(std::remove_if(work.begin(), work.end(), [&](const Particle& part) { return inCone(seed, part); }),  //erase particles from further jet clustering
               work.end());
     }
   }

    if ( useExternalSeeds ) {
      // if (debug_) {dbgCout() << "External seed used!" << std::endl;}
      seeds.erase(seeds.begin());
      // if (debug_){ dbgCout() << "N seeds remaining and: " << seeds.size() << std::endl;}

    }


  }
  return jets;
}
