const unsigned int kR16dAu200_bbcnc    = 0x00000008; // Run16dAu200 BBCS narrow central
const unsigned int kR16dAu200_bbn      = 0x00000010; // Run16dAu200 BBCS narrow 
const unsigned int kR16dAu200_ert4x4b  = 0x00000040; // Run16dAu200 ERT
const unsigned int kR16dAu200_mpcA     = 0x00010000; // Run16dAu200 MPCA
const unsigned int kR16dAu200_mpcB     = 0x00020000; // Run16dAu200 MPCB

const unsigned int kR16dAu39_fvtxbbcnc = 0x00100000; // Run16dAu  FVTXS BBCS central
const unsigned int kR16dAu39_fvtxbbcn  = 0x00400000; // Run16dAu  FVTXS BBCS

//HeavyIons
const Int_t kBinsCen = 60;
const Int_t kBinsVtx = 40;
const Int_t kFlatten = 32;
const Float_t kMinBinVtx = -20.0;
const Float_t kMinBinCen = 0.0;
//const unsigned int kTriggerMask = kR16dAu200_bbcnc | kR16dAu200_bbn; // Run16dAu200 MB
const unsigned int kTriggerMask = kR16dAu200_bbcnc | kR16dAu200_bbn | kR16dAu200_ert4x4b; // Run16dAu200 ERT


//pp
//const Int_t kBinsCen = 1;
//const Int_t kBinsVtx = 40;
//const Int_t kFlatten = 32;
//const Float_t kMinBinVtx = -20.0;
//const Float_t kMinBinCen = 0.0;


//Config for EventChecker
const int kBBCADCMin_S = 71; // Au
const int kBBCADCMax_S = 93; // Au
const int kBBCADCMin_N = 5;  // d
const int kBBCADCMax_N = 50; // d
