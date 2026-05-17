#include "TH1.h"
#include "TGraph.h"
#include "TEfficiency.h"
#include "TH2.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGraphErrors.h"
#include "TVectorD.h"
#include "TTimeStamp.h"
#include <fstream>
#include "TMinuit.h"
#include "TString.h"
#include <vector>
#include <string.h>
#include "TLatex.h"
#include "TPaveStats.h"
#include "TDatime.h"
#include "TColor.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TPrincipal.h"
#include "TDecompChol.h"
#include "TEfficiency.h"
enum Mode { kKO, kAbs };
#include <map>
#include <variant>
 
struct LinearParams {
    double slope, intercept;
};
 
struct GammaParams {
    double slope, intercept, floor, exp;
};
 
struct ModelParamsLinear {
    LinearParams hA, hN, INCL, G4;
};
 
struct ModelParamsGamma {
    GammaParams hA, hN, INCL, G4;
};
 
enum class ProbeType { Pion, Nucleon };
 
struct NucleonTargetParams {
    ModelParamsLinear DifMean;
    ModelParamsLinear DifStd;
    ModelParamsGamma  SumGamma;
};
 
struct PionTargetParams {
    ModelParamsLinear DifMean;
    ModelParamsLinear DifStd;
    ModelParamsLinear SumMean;
    ModelParamsLinear SumStd;
};
 
struct ParticleTargetParams {
    ProbeType probeType;
    std::variant<PionTargetParams, NucleonTargetParams> params;
};
 
// Keyed by [PDG code][target PDG code]
std::map<int, std::map<int, ParticleTargetParams>> fitParams;

void makeAllParams() {

    // PDG 2112 (neutron)
    fitParams[2112][1000060120] = {  // C12
        ProbeType::Nucleon,
        NucleonTargetParams{
            // DifMean
            {
                {-0.0480, -0.4030},  // hA
                {-0.0030, -0.4805},  // hN
                {-0.0991, -0.1968},  // INCL
                {0.0678, -0.4635}  // G4
            },
            // DifStd
            {
                {0.0238, 2.0497},  // hA
                {-0.0146, 1.5578},  // hN
                {-0.0444, 1.0927},  // INCL
                {-0.0296, 1.1317}  // G4
            },
            // SumGamma
            {
                {-0.5662, 0.0000, 0.2119, 2.7880},  // hA
                {-1.4010, 1.5131, 0.0000, 0.1000},  // hN
                {-3.9810, 3.2275, 0.4697, 0.6677},  // INCL
                {-4.2979, 1.1055, 0.3746, 1.3954}  // G4
            }
        }
    };
    fitParams[2112][1000080160] = {  // O16
        ProbeType::Nucleon,
        NucleonTargetParams{
            // DifMean
            {
                {-0.0363, -0.4040},  // hA
                {0.0146, -0.5019},  // hN
                {-0.1060, -0.0903},  // INCL
                {0.0478, -0.3310}  // G4
            },
            // DifStd
            {
                {-0.0022, 2.1431},  // hA
                {-0.0077, 1.6846},  // hN
                {-0.0009, 1.1707},  // INCL
                {-0.0086, 1.1943}  // G4
            },
            // SumGamma
            {
                {-0.0000, 0.1777, 0.0032, 0.9746},  // hA
                {-5.6527, 0.1625, 0.2941, 2.8461},  // hN
                {-4.6053, 9.9999, 0.4538, 0.4330},  // INCL
                {-4.3938, 1.0356, 0.3111, 1.5890}  // G4
            }
        }
    };
    fitParams[2112][1000180400] = {  // Ar40
        ProbeType::Nucleon,
        NucleonTargetParams{
            // DifMean
            {
                {-0.0983, -1.3899},  // hA
                {-0.1625, -0.9717},  // hN
                {-0.3408, -0.8803},  // INCL
                {-0.3147, -1.2028}  // G4
            },
            // DifStd
            {
                {-0.0182, 2.3257},  // hA
                {0.1724, 2.0530},  // hN
                {0.2088, 1.4219},  // INCL
                {0.1920, 1.4619}  // G4
            },
            // SumGamma
            {
                {-0.0239, 0.1237, 0.0006, 4.1641},  // hA
                {-3.8305, 0.2127, 0.1501, 2.0211},  // hN
                {-3.9471, 2.0682, 0.2141, 0.7662},  // INCL
                {-3.2923, 1.8328, 0.0878, 0.7770}  // G4
            }
        }
    };
    fitParams[2112][1000260560] = {  // Fe56
        ProbeType::Nucleon,
        NucleonTargetParams{
            // DifMean
            {
                {-0.0860, -1.0983},  // hA
                {-0.1745, -0.8599},  // hN
                {-0.4010, -0.0581},  // INCL
                {-0.3393, -0.1383}  // G4
            },
            // DifStd
            {
                {-0.0014, 2.4809},  // hA
                {0.2740, 2.1740},  // hN
                {0.2935, 1.4031},  // INCL
                {0.2336, 1.5123}  // G4
            },
            // SumGamma
            {
                {-0.2559, 0.0498, 0.0597, 3.5358},  // hA
                {-4.4328, 0.2101, 0.1197, 2.0547},  // hN
                {-4.1442, 2.9649, 0.1600, 0.6114},  // INCL
                {-3.3238, 3.2311, 0.0000, 0.4775}  // G4
            }
        }
    };

    // PDG 2212 (proton)
    fitParams[2212][1000060120] = {  // C12
        ProbeType::Nucleon,
        NucleonTargetParams{
            // DifMean
            {
                {0.0397, 1.5895},  // hA
                {0.0067, 1.4666},  // hN
                {0.0686, 1.3783},  // INCL
                {-0.1096, 1.7054}  // G4
            },
            // DifStd
            {
                {0.0256, 2.0761},  // hA
                {-0.0241, 1.5651},  // hN
                {-0.0429, 1.0852},  // INCL
                {-0.0383, 1.1547}  // G4
            },
            // SumGamma
            {
                {-0.0000, 0.2097, 0.0031, 4.9989},  // hA
                {-1.5993, 1.7594, 0.0000, 0.1032},  // hN
                {-3.8458, 2.7042, 0.4699, 0.7212},  // INCL
                {-4.4785, 0.9803, 0.3779, 1.6314}  // G4
            }
        }
    };
    fitParams[2212][1000080160] = {  // O16
        ProbeType::Nucleon,
        NucleonTargetParams{
            // DifMean
            {
                {0.0187, 1.5989},  // hA
                {-0.0029, 1.4802},  // hN
                {0.0681, 1.3906},  // INCL
                {-0.1523, 1.8342}  // G4
            },
            // DifStd
            {
                {0.0011, 2.1667},  // hA
                {-0.0089, 1.6803},  // hN
                {0.0092, 1.1461},  // INCL
                {-0.0116, 1.2014}  // G4
            },
            // SumGamma
            {
                {-0.0000, 0.1696, 0.0115, 0.9472},  // hA
                {-4.9704, 7.9247, 0.2394, 0.1825},  // hN
                {-3.5933, 2.4412, 0.4722, 0.8241},  // INCL
                {-4.2475, 1.0173, 0.3022, 1.6352}  // G4
            }
        }
    };
    fitParams[2212][1000180400] = {  // Ar40
        ProbeType::Nucleon,
        NucleonTargetParams{
            // DifMean
            {
                {-0.0086, 0.4249},  // hA
                {-0.1964, 0.9925},  // hN
                {-0.0881, 0.3557},  // INCL
                {-0.5101, 0.8342}  // G4
            },
            // DifStd
            {
                {-0.0230, 2.4054},  // hA
                {0.2043, 2.0487},  // hN
                {0.2477, 1.3776},  // INCL
                {0.1220, 1.5599}  // G4
            },
            // SumGamma
            {
                {-0.0447, 0.1173, 0.0070, 3.2000},  // hA
                {-4.7064, 0.1913, 0.1374, 2.7214},  // hN
                {-4.7041, 6.8098, 0.2036, 0.4407},  // INCL
                {-3.3876, 1.9572, 0.0907, 0.7590}  // G4
            }
        }
    };
    fitParams[2212][1000260560] = {  // Fe56
        ProbeType::Nucleon,
        NucleonTargetParams{
            // DifMean
            {
                {0.0005, 0.7074},  // hA
                {-0.2305, 1.1033},  // hN
                {-0.1581, 1.1520},  // INCL
                {-0.5472, 1.8684}  // G4
            },
            // DifStd
            {
                {-0.0029, 2.5204},  // hA
                {0.3196, 2.1606},  // hN
                {0.3124, 1.3731},  // INCL
                {0.2439, 1.4853}  // G4
            },
            // SumGamma
            {
                {-0.3339, 0.0528, 0.0595, 2.9499},  // hA
                {-4.6425, 0.1992, 0.1069, 2.4347},  // hN
                {-4.0675, 2.9772, 0.1584, 0.5802},  // INCL
                {-3.3271, 3.2049, 0.0000, 0.4905}  // G4
            }
        }
    };

    // PDG -211 (piminus)
    fitParams[-211][1000060120] = {  // C12
        ProbeType::Pion,
        PionTargetParams{
            // DifMean
            {
                {2.7609, -1.8365},  // hA
                {0.1833, -1.3554},  // hN
                {-0.1699, -1.1109},  // INCL
                {0.0545, -0.9973}  // G4
            },
            // DifStd
            {
                {0.2195, 3.4732},  // hA
                {-0.0569, 1.3899},  // hN
                {-0.0584, 0.9972},  // INCL
                {0.1590, 0.9629}  // G4
            },
            // SumMean
            {
                {0.6658, 4.2362},  // hA
                {2.3593, 4.8110},  // hN
                {1.4803, 3.2562},  // INCL
                {2.2654, 1.0971}  // G4
            },
            // SumStd
            {
                {0.6186, 0.8173},  // hA
                {0.3440, 2.7930},  // hN
                {0.8436, 1.4135},  // INCL
                {1.4320, 2.7552}  // G4
            }
        }
    };
    fitParams[-211][1000080160] = {  // O16
        ProbeType::Pion,
        PionTargetParams{
            // DifMean
            {
                {2.7501, -1.9629},  // hA
                {0.0008, -1.1685},  // hN
                {-0.2217, -0.9145},  // INCL
                {0.0345, -0.9206}  // G4
            },
            // DifStd
            {
                {0.2652, 3.4930},  // hA
                {-0.0121, 1.5923},  // hN
                {-0.0230, 1.1189},  // INCL
                {0.2313, 1.0197}  // G4
            },
            // SumMean
            {
                {1.0830, 4.3215},  // hA
                {5.3797, 3.9193},  // hN
                {1.9063, 2.8649},  // INCL
                {4.6220, 1.2158}  // G4
            },
            // SumStd
            {
                {0.9601, 0.9015},  // hA
                {0.7569, 2.9448},  // hN
                {0.9527, 1.3764},  // INCL
                {-0.0833, 3.5171}  // G4
            }
        }
    };
    fitParams[-211][1000180400] = {  // Ar40
        ProbeType::Pion,
        PionTargetParams{
            // DifMean
            {
                {1.6508, -1.8136},  // hA
                {-0.2982, -1.8298},  // hN
                {-0.7577, -1.6395},  // INCL
                {-1.3603, -1.3789}  // G4
            },
            // DifStd
            {
                {0.4772, 3.1882},  // hA
                {0.4976, 2.1689},  // hN
                {0.5777, 1.3554},  // INCL
                {0.2633, 1.6002}  // G4
            },
            // SumMean
            {
                {3.4979, 5.0356},  // hA
                {7.3865, 6.2422},  // hN
                {5.5208, 1.8667},  // INCL
                {8.2553, 2.4079}  // G4
            },
            // SumStd
            {
                {2.3536, 1.5691},  // hA
                {1.2768, 4.3182},  // hN
                {2.9798, 2.6283},  // INCL
                {1.4436, 3.0073}  // G4
            }
        }
    };
    fitParams[-211][1000260560] = {  // Fe56
        ProbeType::Pion,
        PionTargetParams{
            // DifMean
            {
                {1.1089, -1.9992},  // hA
                {-0.2659, -1.6327},  // hN
                {-0.9133, -0.7565},  // INCL
                {-1.1672, -0.3635}  // G4
            },
            // DifStd
            {
                {0.5256, 3.2859},  // hA
                {0.6522, 2.3551},  // hN
                {0.7073, 1.3455},  // INCL
                {0.4936, 1.4834}  // G4
            },
            // SumMean
            {
                {4.8382, 5.8804},  // hA
                {8.8250, 6.4708},  // hN
                {5.9102, 3.4901},  // INCL
                {10.6341, 1.3621}  // G4
            },
            // SumStd
            {
                {2.5376, 2.2535},  // hA
                {1.2986, 4.4985},  // hN
                {1.5502, 2.2818},  // INCL
                {1.7488, 2.9591}  // G4
            }
        }
    };

    // PDG 111 (pizero)
    fitParams[111][1000060120] = {  // C12
        ProbeType::Pion,
        PionTargetParams{
            // DifMean
            {
                {2.9052, -0.1317},  // hA
                {0.0309, 0.4817},  // hN
                {0.0591, 0.5153},  // INCL
                {-0.0937, 0.5965}  // G4
            },
            // DifStd
            {
                {-0.1003, 3.7673},  // hA
                {0.0066, 1.3559},  // hN
                {-0.1664, 1.0968},  // INCL
                {-0.1716, 1.3215}  // G4
            },
            // SumMean
            {
                {0.8463, 4.1349},  // hA
                {1.7051, 5.4571},  // hN
                {1.4390, 3.3941},  // INCL
                {0.4335, 1.9058}  // G4
            },
            // SumStd
            {
                {0.6222, 0.8041},  // hA
                {-0.0207, 3.1001},  // hN
                {1.0085, 1.2471},  // INCL
                {2.2667, 2.0590}  // G4
            }
        }
    };
    fitParams[111][1000080160] = {  // O16
        ProbeType::Pion,
        PionTargetParams{
            // DifMean
            {
                {3.1404, -0.3688},  // hA
                {0.0321, 0.4539},  // hN
                {-0.0321, 0.6271},  // INCL
                {-0.1328, 0.7142}  // G4
            },
            // DifStd
            {
                {0.1547, 3.6580},  // hA
                {-0.0448, 1.6472},  // hN
                {-0.0366, 1.2038},  // INCL
                {-0.2173, 1.4185}  // G4
            },
            // SumMean
            {
                {1.3092, 4.1830},  // hA
                {3.1107, 5.4940},  // hN
                {1.7053, 3.0618},  // INCL
                {3.2341, 1.9929}  // G4
            },
            // SumStd
            {
                {0.9129, 0.9172},  // hA
                {0.6364, 2.9601},  // hN
                {0.8161, 1.4066},  // INCL
                {0.4365, 3.1033}  // G4
            }
        }
    };
    fitParams[111][1000180400] = {  // Ar40
        ProbeType::Pion,
        PionTargetParams{
            // DifMean
            {
                {2.0891, -0.5501},  // hA
                {-0.6140, -0.0622},  // hN
                {-0.5836, -0.3096},  // INCL
                {-1.1870, 0.0335}  // G4
            },
            // DifStd
            {
                {0.4758, 3.2005},  // hA
                {0.4925, 2.1847},  // hN
                {0.5297, 1.4072},  // INCL
                {0.0295, 1.7338}  // G4
            },
            // SumMean
            {
                {3.4947, 4.9809},  // hA
                {6.8813, 6.9028},  // hN
                {5.6032, 1.6268},  // INCL
                {9.2720, 0.8899}  // G4
            },
            // SumStd
            {
                {2.2506, 1.6529},  // hA
                {1.6484, 4.0574},  // hN
                {3.2671, 2.6175},  // INCL
                {1.3174, 3.4054}  // G4
            }
        }
    };
    fitParams[111][1000260560] = {  // Fe56
        ProbeType::Pion,
        PionTargetParams{
            // DifMean
            {
                {1.4231, -0.5136},  // hA
                {-0.5220, -0.0200},  // hN
                {-0.6510, 0.5007},  // INCL
                {-0.9026, 0.9184}  // G4
            },
            // DifStd
            {
                {0.3810, 3.4108},  // hA
                {0.5532, 2.4634},  // hN
                {0.7027, 1.3688},  // INCL
                {0.2354, 1.6588}  // G4
            },
            // SumMean
            {
                {4.9143, 5.6408},  // hA
                {8.0224, 7.3538},  // hN
                {5.6536, 3.6493},  // INCL
                {10.8326, 0.7434}  // G4
            },
            // SumStd
            {
                {2.6421, 2.2420},  // hA
                {1.7590, 4.3464},  // hN
                {1.6016, 2.1765},  // INCL
                {1.8533, 3.0423}  // G4
            }
        }
    };

    // PDG 211 (piplus)
    fitParams[211][1000060120] = {  // C12
        ProbeType::Pion,
        PionTargetParams{
            // DifMean
            {
                {2.9002, 1.7511},  // hA
                {-0.1257, 2.3318},  // hN
                {0.1201, 2.2436},  // INCL
                {-0.2030, 2.1775}  // G4
            },
            // DifStd
            {
                {-0.1002, 3.7254},  // hA
                {-0.0056, 1.3204},  // hN
                {-0.1042, 1.0109},  // INCL
                {0.1389, 0.9416}  // G4
            },
            // SumMean
            {
                {1.0770, 4.0994},  // hA
                {2.1249, 4.8008},  // hN
                {1.3061, 3.4174},  // INCL
                {0.7052, 2.0391}  // G4
            },
            // SumStd
            {
                {0.6136, 0.8361},  // hA
                {0.2806, 3.2332},  // hN
                {0.8760, 1.3757},  // INCL
                {1.2589, 2.7659}  // G4
            }
        }
    };
    fitParams[211][1000080160] = {  // O16
        ProbeType::Pion,
        PionTargetParams{
            // DifMean
            {
                {2.8132, 1.9100},  // hA
                {-0.0463, 2.1889},  // hN
                {0.1818, 2.1437},  // INCL
                {-0.1161, 2.2395}  // G4
            },
            // DifStd
            {
                {-0.1951, 3.9639},  // hA
                {0.0106, 1.5698},  // hN
                {0.0372, 1.0794},  // INCL
                {0.1149, 1.0644}  // G4
            },
            // SumMean
            {
                {1.4981, 4.2037},  // hA
                {3.2590, 5.0669},  // hN
                {1.8017, 2.9937},  // INCL
                {2.5676, 3.1655}  // G4
            },
            // SumStd
            {
                {0.9073, 0.9469},  // hA
                {0.7763, 3.4013},  // hN
                {0.8655, 1.4225},  // INCL
                {1.0681, 2.2491}  // G4
            }
        }
    };
    fitParams[211][1000180400] = {  // Ar40
        ProbeType::Pion,
        PionTargetParams{
            // DifMean
            {
                {2.2903, 0.8905},  // hA
                {-0.6633, 1.4800},  // hN
                {-0.4752, 1.1072},  // INCL
                {-1.0096, 1.3981}  // G4
            },
            // DifStd
            {
                {0.2987, 3.2604},  // hA
                {0.3182, 2.2811},  // hN
                {0.7081, 1.2767},  // INCL
                {0.2957, 1.3652}  // G4
            },
            // SumMean
            {
                {3.7611, 4.9256},  // hA
                {6.6766, 7.1275},  // hN
                {4.9647, 1.8701},  // INCL
                {8.6358, 1.7998}  // G4
            },
            // SumStd
            {
                {2.1374, 1.6616},  // hA
                {1.2266, 4.5166},  // hN
                {3.7494, 2.4962},  // INCL
                {1.4930, 3.0141}  // G4
            }
        }
    };
    fitParams[211][1000260560] = {  // Fe56
        ProbeType::Pion,
        PionTargetParams{
            // DifMean
            {
                {1.8171, 0.9453},  // hA
                {-0.7655, 1.6732},  // hN
                {-0.6367, 1.9307},  // INCL
                {-0.6924, 2.2437}  // G4
            },
            // DifStd
            {
                {0.3591, 3.3994},  // hA
                {0.6809, 2.3669},  // hN
                {0.7452, 1.2930},  // INCL
                {0.3858, 1.4334}  // G4
            },
            // SumMean
            {
                {5.0295, 5.7037},  // hA
                {9.1261, 6.7006},  // hN
                {5.6704, 3.6514},  // INCL
                {10.2044, 1.6007}  // G4
            },
            // SumStd
            {
                {2.5278, 2.2597},  // hA
                {1.7786, 4.6312},  // hN
                {1.6104, 2.1473},  // INCL
                {2.0127, 2.6716}  // G4
            }
        }
    };
}
const auto& getPion(int pdg, int target) {
    return std::get<PionTargetParams>(fitParams.at(pdg).at(target).params);
}

const auto& getNucleon(int pdg, int target) {
    return std::get<NucleonTargetParams>(fitParams.at(pdg).at(target).params);
}
double GetVisEReweight(
    TH2D* hist_nom,
    TH2D* hist_alt,
    double KEini,
    double Ebias)
{
    if (!hist_nom || !hist_alt) return 1.;

    int idx_KEini = hist_nom->GetXaxis()->FindBin(KEini);
    int idx_Ebias = hist_nom->GetYaxis()->FindBin(Ebias);

    double weight_nom = hist_nom->GetBinContent(idx_KEini, idx_Ebias);
    double weight_alt = hist_alt->GetBinContent(idx_KEini, idx_Ebias);

    if (weight_nom == 0.) return 1.;


    return weight_alt / weight_nom;
}
TH1D* makeHistogram(TTree* t, int pdg, std::string label, Mode mode, bool useThreshold,int A1, int target, bool compound, bool RW, std::string pstr="pionp",std::string reweightModel="INCL")
{

    int min=0; int max=40;
    if (A1<40){ min=0; max=A1+2;}
    TH1D *hist = new TH1D(Form("hist_%s",label.c_str()),
                          Form("hist_%s",label.c_str()),max,min,max);
    if (!t || t->GetEntries() < 1) return hist;

    int pdgh[100], nh; 
    double Eh[100], mh[100];
    int probe_fsi;
    int npi0, npip, npim;
    double KEini, Eini;
    
    t->SetBranchAddress("pdgh",      pdgh);
    t->SetBranchAddress("npi0",   &npi0);
    t->SetBranchAddress("npip",   &npip);
    t->SetBranchAddress("npim",   &npim);
    t->SetBranchAddress("ke",       &KEini);
    t->SetBranchAddress("e",        &Eini);
    t->SetBranchAddress("nh",       &nh);
    t->SetBranchAddress("Eh",        Eh);
    t->SetBranchAddress("mh",        mh);
    t->SetBranchAddress("probe_fsi",&probe_fsi);

    
    int required_fsi = (mode == kAbs) ? 5 : 6;
    double thresh = useThreshold ? 0.005 : 0.0;

    for (Long64_t j=0; j<t->GetEntries(); j++){
        t->GetEntry(j);
       int mult=0;
       int piMult=0;
       int diff=0;
    double visE=0;

       for(int i=0; i<nh; i++){
            int apdg = abs(pdgh[i]);
            double KE = Eh[i] - mh[i];
            if (apdg==211|| apdg==111){
                piMult++;
            }
            
            if (KE <= thresh) continue;

            if (apdg==2212 || apdg==2112){
                mult++;
                if (apdg==2212) diff++;
                else diff--;
            }
           
            if (apdg > 1000000000){
                int Z = (apdg/10000) % 1000;
                int A = (apdg/10) % 1000;
                int N = A - Z;

                   if (Z > 2) continue;
             if (KE<=thresh*A) continue;

             if (compound==true)    mult += Z + N;
            if (compound==true) diff+=Z-N;
            }
            
        }
        if (required_fsi==6 && (mult<3 || piMult!=0)) continue;
        if (required_fsi==5 && (piMult!=0 || mult<2)) continue;


        for (int i = 0; i < nh; i++) {           // fixed: nf/pdgf/Ef -> nh/pdgh/Eh
            int    apdg = abs(pdgh[i]);
            double KE   = Eh[i] - mh[i];

            if (apdg == 211) {
                visE += Eh[i];                    // fixed: Ef -> Eh
            }
            if (apdg == 22) {
                //if (KE <= 0.5 * thresh) continue;
                visE += Eh[i];                    // fixed: Ef -> Eh

            }
            
            if (apdg == 111) {
                visE += Eh[i];                    // fixed: Ef -> Eh


            }
            if (KE <= thresh) continue;

            if (apdg == 2212 || apdg == 2112) {
                if (apdg == 2212) visE += KE;
            }

            if (apdg > 1000000000) {
                int Z = (apdg / 10000) % 1000;
                int A = (apdg / 10)    % 1000;
                int N = A - Z;

                if (Z > 2)          continue;
                if (KE <= thresh*A) continue;


                   if (compound==1){
                                        visE += KE;}
            }
        }

        double Ehad = visE;                       // fixed: removed re-declaration of visE, named clearly
        double Ebias = 0.0;
        if (mode==kKO)
            Ebias = (KEini - Ehad) / KEini;
        else 
            Ebias = (Eini - Ehad) / Eini;

if (RW) {
    TFile* f3=TFile::Open("FSI_KOAbs_Evis_reweight_template.root","READ");
    TFile* f2 = TFile::Open("FSI_KOAbs_hA2018rwNucl_templates.root", "READ");
        TH2D* hist_nom = (TH2D*)f2->Get(Form("rw%s_%s_KEini_vs_Ebias",reweightModel.c_str(), pstr.c_str()));
    TH2D* hist_alt = (TH2D*)f3->Get(Form("%s_%s_KEini_vs_Ebias",     reweightModel.c_str(), pstr.c_str()));
    double visERW=GetVisEReweight(hist_nom, hist_alt, KEini, Ebias);
    double hAEst=1;
    double AltEst=1;
    if (mode==kAbs){
        const auto& tp = getPion(pdg, target);


        const LinearParams& altSumMean = (reweightModel=="hN2018") ? tp.SumMean.hN   :
                                         (reweightModel=="Geant4")  ? tp.SumMean.G4   :
                                                                       tp.SumMean.INCL;
        const LinearParams& altSumStd  = (reweightModel=="hN2018") ? tp.SumStd.hN    :
                                         (reweightModel=="Geant4")  ? tp.SumStd.G4    :
                                                                       tp.SumStd.INCL;

        double meanAltSum = altSumMean.slope * KEini + altSumMean.intercept;
        double stdAltSum  = altSumStd.slope  * KEini + altSumStd.intercept;
        double meanhASum  = tp.SumMean.hA.slope * KEini + tp.SumMean.hA.intercept;
        double stdhASum   = tp.SumStd.hA.slope  * KEini + tp.SumStd.hA.intercept;



        const LinearParams& altDifMean = (reweightModel=="hN2018") ? tp.DifMean.hN   :
                                         (reweightModel=="Geant4")  ? tp.DifMean.G4   :
                                                                       tp.DifMean.INCL;
        const LinearParams& altDifStd  = (reweightModel=="hN2018") ? tp.DifStd.hN    :
                                         (reweightModel=="Geant4")  ? tp.DifStd.G4    :
                                                                       tp.DifStd.INCL;

        double meanAltDif = altDifMean.slope * KEini + altDifMean.intercept;
        double stdAltDif  = altDifStd.slope  * KEini + altDifStd.intercept;
        double meanhADif  = tp.DifMean.hA.slope * KEini + tp.DifMean.hA.intercept;
        double stdhADif   = tp.DifStd.hA.slope  * KEini + tp.DifStd.hA.intercept;

        double hAEstSum  = TMath::Gaus(mult, meanhASum,  stdhASum,  1);
        double AltEstSum = TMath::Gaus(mult, meanAltSum, stdAltSum, 1);
        double hAEstDif  = TMath::Gaus(diff, meanhADif,  stdhADif,  1);
        double AltEstDif = TMath::Gaus(diff, meanAltDif, stdAltDif, 1);
        hAEst=hAEstSum*hAEstDif; AltEst=AltEstSum*AltEstDif;
    
    }
    else{


        const auto& tp = getNucleon(pdg, target);

        const GammaParams& altGamma = (reweightModel=="hN2018") ? tp.SumGamma.hN   :
                                      (reweightModel=="Geant4")  ? tp.SumGamma.G4   :
                                                                   tp.SumGamma.INCL;

        double gammaAlt = altGamma.intercept * TMath::Exp(TMath::Power(KEini, altGamma.exp) * altGamma.slope) + altGamma.floor;
        double gammaHA  = tp.SumGamma.hA.intercept * TMath::Exp(TMath::Power(KEini, tp.SumGamma.hA.exp) * tp.SumGamma.hA.slope) + tp.SumGamma.hA.floor;
        hAEst  = TMath::Exp(-(mult - 3) * gammaHA);
        AltEst = TMath::Exp(-(mult - 3) * gammaAlt);
    }
    double rw     = visERW*(AltEst / hAEst);
    if(rw<0.001) rw=0.001; if (rw>100) rw=100;
           
    
        hist->Fill(mult, rw);
        } else {
            hist->Fill(mult, 1);
        }
    }

    return hist;
}

// =======================================================
// Main plotting function
// =======================================================
void plotTotNuclwAllRW(int probe=211, double ke=0.175, int target=1000180400,
                  bool useThreshold = true, std::string reweightModel="hN2018"){
makeAllParams();
    std::string probeStr;
    Mode mode = (probe==2212 || probe==-2212 || probe==2112) ? kKO : kAbs;
    bool isPion = (probe==211 || probe==-211 || probe==111);
    // ---------------------------------------------------
    // Probe string (for filenames)
    // ---------------------------------------------------
    if (mode == kAbs){
        if (probe==211) probeStr="piPlus";
        if (probe==-211) probeStr="piMinus";
        if (probe==111) probeStr="pi0";
    } else {
        if (probe==2212) probeStr="protonPlus";
        if (probe==-2212) probeStr="protonMinus";
        if (probe==2112) probeStr="neutron";
    }

    std::string keStr = std::to_string(ke).substr(0,5);
    std::string keLow=std::to_string(ke-0.025).substr(0, 4);
    std::string keHigh=std::to_string(ke+0.025).substr(0, 4);
    gROOT->LoadMacro("protoDUNEStyle.C");
    gROOT->SetStyle("protoDUNEStyle");
    gROOT->ForceStyle();
    gROOT->ForceStyle();
    gStyle->SetTitleX(0.45);
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);

  TCanvas *c1 = new TCanvas("c1","c1",200,10,700,500);

    // ---------------------------------------------------
    // Load files
    // ---------------------------------------------------
    TFile *f1 = new TFile(Form("$PSCRATCH/hadronFiles/%s_%d_%sGeV_hA2018_100k.ginuke.root",
                              probeStr.c_str(),target,keStr.c_str()));
    TFile *f2 = new TFile(Form("$PSCRATCH/hadronFiles/%s_%d_%sGeV_hN2018_100k.ginuke.root",
                              probeStr.c_str(),target,keStr.c_str()));
    TFile *f3 = new TFile(Form("$PSCRATCH/hadronFiles/%s_%d_%sGeV_HINCL_100k.ginuke.root",
                              probeStr.c_str(),target,keStr.c_str()));
    TFile *f4 = new TFile(Form("$PSCRATCH/hadronFiles/%s_%d_%sGeV_HG4BertCasc_100k.ginuke.root",probeStr.c_str(),target,keStr.c_str()));

    int Z = (target/10000) % 1000;
    int A = (target/10) % 1000;
    double min=-20; double max=20;
    if (A<40){ min=0; max=12;}
                              

    TTree *t1 = (TTree*)f1->Get("ginuke");
    TTree *t2 = (TTree*)f2->Get("ginuke");
    TTree *t3 = (TTree*)f3->Get("ginuke");
    TTree *t4 = (TTree*)f4->Get("ginuke");

    // ---------------------------------------------------
    // Histograms
    // ---------------------------------------------------
    int compound=0;
    TH1D* h1 = makeHistogram(t1,probe,"20i",mode,useThreshold,A,target,compound,1,reweightModel);
    TH1D* h2 = makeHistogram(t2,probe,"20j",mode,useThreshold,A,target,compound,0);
    TH1D* h3 = makeHistogram(t3,probe,"20k",mode,useThreshold,A,target,compound,0);
    TH1D* h4 = makeHistogram(t4,probe,"20l",mode,useThreshold,A,target,compound,0);

h1->Scale(1.0 / h1->Integral());
h2->Scale(1.0 / h2->Integral());
h3->Scale(1.0 / h3->Integral());
h4->Scale(1.0 / h4->Integral());

    h1->SetLineColor(kBlack);
    h2->SetLineColor(kRed);   h2->SetLineStyle(2);
    h3->SetLineColor(kBlue);  h3->SetLineStyle(4);
    h4->SetLineColor(kGreen); h4->SetLineStyle(3);

    // ---------------------------------------------------
    // Labels and titles
    // ---------------------------------------------------
    std::string title;
    if (mode==kAbs)
        h1->SetTitle("Pion Absorption");
    else
        h1->SetTitle("Nucleon Knockout");


    if (useThreshold)
        h1->GetXaxis()->SetTitle("N_{p}+N_{n} (T_{p,n}>5 MeV)");
    else
        h1->GetXaxis()->SetTitle("N_{p}+N_{n}");
    
    h1->GetYaxis()->SetTitle("Fraction of Absorption Int.");
    if (mode!=kAbs)    h1->GetYaxis()->SetTitle("Fraction of Knockout Int.");
    h1->GetXaxis()->CenterTitle();
    h1->GetYaxis()->CenterTitle();
    if (A==12) h1->GetYaxis()->SetRangeUser(0,1);
    else h1->GetYaxis()->SetRangeUser(0,0.6);
    if (h4->GetMaximum()>0.5) h1->GetYaxis()->SetRangeUser(0,1);
    h1->GetYaxis()->SetTitleOffset(1.3);

    h1->Draw("HIST");
    h2->Draw("HIST SAME");
    h3->Draw("HIST SAME");
    h4->Draw("HIST SAME");

    // ---------------------------------------------------
    // Legend header (correct particle + target)
    // ---------------------------------------------------
    std::string particleLabel;

    if (mode == kAbs){
        if (probe==211) particleLabel = "#pi^{+}";
        if (probe==-211) particleLabel = "#pi^{-}";
        if (probe==111) particleLabel = "#pi^{0}";
    } else {
        if (probe==2212) particleLabel = "p^{+}";
        if (probe==-2212) particleLabel = "p^{-}";
        if (probe==2112) particleLabel = "n";
    }

    std::string targetStr = "C";
    if (target==1000822080) targetStr="Pb";
    if (target==1000260560) targetStr="Fe";
    if (target==1000080160) targetStr="O";
    if (target==1000180400) targetStr="Ar";


    TLegend *leg = new TLegend(0.55, 0.65, 0.45, 0.88);
leg->SetHeader(Form("%s+%s (%s#font[122]{-}%s GeV)",
                    particleLabel.c_str(),
                    targetStr.c_str(),
                    keLow.c_str(), keHigh.c_str()));
leg->SetTextFont(42);       // standard scalable font
leg->SetTextSize(0.05);     // smaller, fits entries



leg->AddEntry(h1,Form("hA2018 Reweighted to %s",reweightModel.c_str()),"l");
leg->AddEntry(h2,"hN2018","l");
leg->AddEntry(h3,"INCL++","l");
leg->AddEntry(h4,"Geant4","l");

leg->Draw();

    c1->SetLeftMargin(0.13);
    c1->SetRightMargin(0.025);
    c1->Update();  // commit before drawing

    // ---------------------------------------------------
    // Output naming (NO "merged")
    // ---------------------------------------------------
    std::string channelStr="nucleonSumThreshPionAbs";

    if (mode == kAbs)
        channelStr = useThreshold ?
            "nucleonSumThreshPionAbs" :
            "nucleonSumPionAbs";
    if (mode==kKO)
        channelStr = useThreshold ?
            "nucleonSumThreshKO" :
            "nucleonSumKO";
    // Argon Sum parameters




    c1->Print(Form("../plotting/%s_%s_%s_%d_%d_%dRWAll%s.pdf",
                  channelStr.c_str(),
                  probeStr.c_str(),
                  keStr.c_str(),
                  target,useThreshold,compound,reweightModel.c_str()));



}