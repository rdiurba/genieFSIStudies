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
                {-0.0458, -0.4058},  // hA
                {0.0098, -0.4984},  // hN
                {-0.0605, -0.2798},  // INCL
                {0.0709, -0.4615}  // G4
            },
            // DifStd
            {
                {0.0131, 2.0646},  // hA
                {-0.0168, 1.5610},  // hN
                {-0.0610, 0.9932},  // INCL
                {-0.0460, 1.1354}  // G4
            },
            // SumGamma
            {
                {-0.0004, 0.0176, 0.1985, 0.1000},  // hA
                {-1.3001, 1.3636, 0.0000, 0.1001},  // hN
                {-1.7230, 1.1991, 0.0000, 0.1373},  // INCL
                {-5.9100, 0.3582, 0.2882, 3.0541}  // G4
            }
        }
    };
    fitParams[2112][1000080160] = {  // O16
        ProbeType::Nucleon,
        NucleonTargetParams{
            // DifMean
            {
                {-0.0138, -0.4352},  // hA
                {0.0292, -0.5224},  // hN
                {-0.0759, -0.1609},  // INCL
                {0.0528, -0.3374}  // G4
            },
            // DifStd
            {
                {-0.0045, 2.1463},  // hA
                {-0.0204, 1.7025},  // hN
                {-0.0538, 1.1200},  // INCL
                {-0.0425, 1.2290}  // G4
            },
            // SumGamma
            {
                {-0.0131, 0.0066, 0.1803, 0.1002},  // hA
                {-1.5255, 1.4372, 0.0000, 0.1053},  // hN
                {-7.6522, 9.8971, 0.2289, 0.6249},  // INCL
                {-6.3157, 0.4936, 0.2881, 2.7992}  // G4
            }
        }
    };
    fitParams[2112][1000180400] = {  // Ar40
        ProbeType::Nucleon,
        NucleonTargetParams{
            // DifMean
            {
                {-0.0498, -1.4574},  // hA
                {-0.1262, -1.0222},  // hN
                {-0.2771, -0.9959},  // INCL
                {-0.3411, -1.1834}  // G4
            },
            // DifStd
            {
                {-0.0045, 2.3063},  // hA
                {0.1143, 2.1339},  // hN
                {0.1603, 1.4664},  // INCL
                {0.1468, 1.5241}  // G4
            },
            // SumGamma
            {
                {-1.7379, 0.6495, 0.0000, 0.1003},  // hA
                {-6.2676, 6.9147, 0.1467, 0.4568},  // hN
                {-6.6095, 6.2220, 0.1407, 0.6269},  // INCL
                {-3.3357, 0.8267, 0.0913, 1.3471}  // G4
            }
        }
    };
    fitParams[2112][1000260560] = {  // Fe56
        ProbeType::Nucleon,
        NucleonTargetParams{
            // DifMean
            {
                {-0.0384, -1.1646},  // hA
                {-0.1299, -0.9215},  // hN
                {-0.3301, -0.1915},  // INCL
                {-0.3707, -0.1054}  // G4
            },
            // DifStd
            {
                {0.0084, 2.4672},  // hA
                {0.2016, 2.2740},  // hN
                {0.2428, 1.4377},  // INCL
                {0.1860, 1.5673}  // G4
            },
            // SumGamma
            {
                {-1.7661, 0.5261, 0.0000, 0.1838},  // hA
                {-6.0021, 0.1411, 0.1204, 3.3621},  // hN
                {-4.8448, 0.3332, 0.0985, 2.1657},  // INCL
                {-3.2130, 0.6480, 0.0569, 1.5080}  // G4
            }
        }
    };

    // PDG 2212 (proton)
    fitParams[2212][1000060120] = {  // C12
        ProbeType::Nucleon,
        NucleonTargetParams{
            // DifMean
            {
                {0.0181, 1.6195},  // hA
                {-0.0048, 1.4827},  // hN
                {0.0520, 1.4161},  // INCL
                {-0.0790, 1.6617}  // G4
            },
            // DifStd
            {
                {0.0265, 2.0748},  // hA
                {-0.0193, 1.5583},  // hN
                {-0.0633, 0.9864},  // INCL
                {-0.0515, 1.1469}  // G4
            },
            // SumGamma
            {
                {-0.0021, 0.1255, 0.0917, 0.1041},  // hA
                {-1.5328, 1.6431, 0.0000, 0.1021},  // hN
                {-1.7220, 1.2051, 0.0000, 0.1248},  // INCL
                {-6.1557, 7.3755, 0.2892, 0.6325}  // G4
            }
        }
    };
    fitParams[2212][1000080160] = {  // O16
        ProbeType::Nucleon,
        NucleonTargetParams{
            // DifMean
            {
                {-0.0034, 1.6298},  // hA
                {0.0010, 1.4747},  // hN
                {0.0615, 1.4178},  // INCL
                {-0.1070, 1.7680}  // G4
            },
            // DifStd
            {
                {-0.0059, 2.1764},  // hA
                {-0.0112, 1.6833},  // hN
                {-0.0580, 1.1150},  // INCL
                {-0.0348, 1.2180}  // G4
            },
            // SumGamma
            {
                {-0.0001, 0.0100, 0.1771, 0.1000},  // hA
                {-1.6654, 1.5777, 0.0000, 0.1215},  // hN
                {-6.8995, 4.0652, 0.2352, 0.7567},  // INCL
                {-5.6316, 0.4995, 0.2827, 2.7168}  // G4
            }
        }
    };
    fitParams[2212][1000180400] = {  // Ar40
        ProbeType::Nucleon,
        NucleonTargetParams{
            // DifMean
            {
                {-0.0069, 0.4224},  // hA
                {-0.1789, 0.9682},  // hN
                {-0.1022, 0.3722},  // INCL
                {-0.4763, 0.7756}  // G4
            },
            // DifStd
            {
                {-0.0101, 2.3873},  // hA
                {0.1621, 2.1070},  // hN
                {0.1980, 1.4204},  // INCL
                {0.0903, 1.5971}  // G4
            },
            // SumGamma
            {
                {-1.8203, 0.7022, 0.0000, 0.1014},  // hA
                {-6.1361, 5.6958, 0.1339, 0.5170},  // hN
                {-6.2255, 4.0756, 0.1401, 0.6827},  // INCL
                {-3.4799, 0.9218, 0.0915, 1.2652}  // G4
            }
        }
    };
    fitParams[2212][1000260560] = {  // Fe56
        ProbeType::Nucleon,
        NucleonTargetParams{
            // DifMean
            {
                {-0.0012, 0.7096},  // hA
                {-0.2181, 1.0861},  // hN
                {-0.1576, 1.1421},  // INCL
                {-0.5327, 1.8451}  // G4
            },
            // DifStd
            {
                {0.0088, 2.5039},  // hA
                {0.2508, 2.2552},  // hN
                {0.2552, 1.4138},  // INCL
                {0.1849, 1.5511}  // G4
            },
            // SumGamma
            {
                {-1.7895, 0.5399, 0.0000, 0.1887},  // hA
                {-6.4229, 6.6016, 0.1044, 0.4978},  // hN
                {-4.8679, 0.3490, 0.0993, 2.0361},  // INCL
                {-3.2558, 0.7821, 0.0525, 1.3369}  // G4
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
                {-0.1156, -1.2158},  // INCL
                {0.0579, -1.0079}  // G4
            },
            // DifStd
            {
                {0.2195, 3.4732},  // hA
                {-0.0569, 1.3899},  // hN
                {-0.1274, 0.9037},  // INCL
                {0.1882, 0.9306}  // G4
            },
            // SumMean
            {
                {0.6658, 4.2362},  // hA
                {2.3593, 4.8110},  // hN
                {1.6460, 3.8559},  // INCL
                {0.4414, 2.0446}  // G4
            },
            // SumStd
            {
                {0.6186, 0.8173},  // hA
                {0.3440, 2.7930},  // hN
                {13.1251, -0.9103},  // INCL
                {2.7875, 2.9880}  // G4
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
                {-0.2474, -0.9797},  // INCL
                {0.0322, -0.9268}  // G4
            },
            // DifStd
            {
                {0.2652, 3.4930},  // hA
                {-0.0121, 1.5923},  // hN
                {-0.1428, 1.0599},  // INCL
                {0.2114, 1.0149}  // G4
            },
            // SumMean
            {
                {1.0830, 4.3215},  // hA
                {5.3797, 3.9193},  // hN
                {-0.1347, 2.7288},  // INCL
                {3.9052, 0.3812}  // G4
            },
            // SumStd
            {
                {0.9601, 0.9015},  // hA
                {0.7569, 2.9448},  // hN
                {9.5732, 3.5767},  // INCL
                {3.6217, 2.6746}  // G4
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
                {-0.8457, -1.6877},  // INCL
                {-1.3908, -1.3807}  // G4
            },
            // DifStd
            {
                {0.4772, 3.1882},  // hA
                {0.4976, 2.1689},  // hN
                {0.5290, 1.3862},  // INCL
                {0.2569, 1.6061}  // G4
            },
            // SumMean
            {
                {3.4979, 5.0356},  // hA
                {7.3865, 6.2422},  // hN
                {8.5453, 2.4680},  // INCL
                {11.7655, 0.3049}  // G4
            },
            // SumStd
            {
                {2.3536, 1.5691},  // hA
                {1.2768, 4.3182},  // hN
                {1.8897, 4.6942},  // INCL
                {2.3445, 3.6421}  // G4
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
                {-0.9714, -0.7995},  // INCL
                {-1.2066, -0.3517}  // G4
            },
            // DifStd
            {
                {0.5256, 3.2859},  // hA
                {0.6522, 2.3551},  // hN
                {0.6318, 1.3908},  // INCL
                {0.4715, 1.4936}  // G4
            },
            // SumMean
            {
                {4.8382, 5.8804},  // hA
                {8.8250, 6.4708},  // hN
                {10.2677, 3.3301},  // INCL
                {14.4090, -0.7734}  // G4
            },
            // SumStd
            {
                {2.5376, 2.2535},  // hA
                {1.2986, 4.4985},  // hN
                {2.4414, 4.1590},  // INCL
                {2.9453, 3.4111}  // G4
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
                {0.0060, 0.5442},  // INCL
                {-0.1184, 0.6118}  // G4
            },
            // DifStd
            {
                {-0.1003, 3.7673},  // hA
                {0.0066, 1.3559},  // hN
                {-0.1776, 0.9585},  // INCL
                {-0.1778, 1.3106}  // G4
            },
            // SumMean
            {
                {0.8463, 4.1349},  // hA
                {1.7051, 5.4571},  // hN
                {0.5703, 5.0099},  // INCL
                {0.0003, 2.0002}  // G4
            },
            // SumStd
            {
                {0.6222, 0.8041},  // hA
                {-0.0207, 3.1001},  // hN
                {4.0335, 3.7402},  // INCL
                {3.8872, 2.0507}  // G4
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
                {-0.0558, 0.6370},  // INCL
                {-0.1467, 0.7228}  // G4
            },
            // DifStd
            {
                {0.1547, 3.6580},  // hA
                {-0.0448, 1.6472},  // hN
                {-0.1356, 1.1183},  // INCL
                {-0.2467, 1.4226}  // G4
            },
            // SumMean
            {
                {1.3092, 4.1830},  // hA
                {3.1107, 5.4940},  // hN
                {2.4582, 1.3847},  // INCL
                {0.9516, 2.1067}  // G4
            },
            // SumStd
            {
                {0.9129, 0.9172},  // hA
                {0.6364, 2.9601},  // hN
                {21.0099, -0.2930},  // INCL
                {2.9113, 2.6140}  // G4
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
                {-0.6430, -0.3426},  // INCL
                {-1.2156, 0.0380}  // G4
            },
            // DifStd
            {
                {0.4758, 3.2005},  // hA
                {0.4925, 2.1847},  // hN
                {0.4928, 1.4221},  // INCL
                {0.0108, 1.7389}  // G4
            },
            // SumMean
            {
                {3.4947, 4.9809},  // hA
                {6.8813, 6.9028},  // hN
                {8.4174, 2.7382},  // INCL
                {15.0807, -3.5560}  // G4
            },
            // SumStd
            {
                {2.2506, 1.6529},  // hA
                {1.6484, 4.0574},  // hN
                {1.7047, 4.7834},  // INCL
                {2.4725, 4.1813}  // G4
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
                {-0.7077, 0.4933},  // INCL
                {-0.9195, 0.9232}  // G4
            },
            // DifStd
            {
                {0.3810, 3.4108},  // hA
                {0.5532, 2.4634},  // hN
                {0.6658, 1.3784},  // INCL
                {0.2116, 1.6657}  // G4
            },
            // SumMean
            {
                {4.9143, 5.6408},  // hA
                {8.0224, 7.3538},  // hN
                {10.7279, 3.1948},  // INCL
                {15.1094, -2.1593}  // G4
            },
            // SumStd
            {
                {2.6421, 2.2420},  // hA
                {1.7590, 4.3464},  // hN
                {2.9992, 3.8938},  // INCL
                {3.3541, 3.4803}  // G4
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
                {0.0878, 2.3290},  // INCL
                {-0.2527, 2.2158}  // G4
            },
            // DifStd
            {
                {-0.1002, 3.7254},  // hA
                {-0.0056, 1.3204},  // hN
                {-0.1019, 0.8634},  // INCL
                {0.1635, 0.9133}  // G4
            },
            // SumMean
            {
                {1.0699, 4.1030},  // hA
                {2.0520, 4.8404},  // hN
                {-1.9430, 6.1858},  // INCL
                {0.2475, 1.9571}  // G4
            },
            // SumStd
            {
                {0.6212, 0.8323},  // hA
                {0.3243, 3.2116},  // hN
                {6.9417, 2.1470},  // INCL
                {2.1653, 3.2065}  // G4
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
                {0.1690, 2.2358},  // INCL
                {-0.1429, 2.2621}  // G4
            },
            // DifStd
            {
                {-0.1951, 3.9639},  // hA
                {0.0106, 1.5698},  // hN
                {-0.1041, 1.0261},  // INCL
                {0.1065, 1.0518}  // G4
            },
            // SumMean
            {
                {1.4981, 4.2037},  // hA
                {3.2590, 5.0669},  // hN
                {1.7280, 1.4570},  // INCL
                {1.6754, 1.6791}  // G4
            },
            // SumStd
            {
                {0.9073, 0.9469},  // hA
                {0.7763, 3.4013},  // hN
                {13.2244, 2.3207},  // INCL
                {1.5898, 3.6865}  // G4
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
                {-0.4911, 1.0835},  // INCL
                {-1.0216, 1.3981}  // G4
            },
            // DifStd
            {
                {0.2987, 3.2604},  // hA
                {0.3182, 2.2811},  // hN
                {0.6505, 1.2921},  // INCL
                {0.2593, 1.3767}  // G4
            },
            // SumMean
            {
                {3.7611, 4.9256},  // hA
                {6.6766, 7.1275},  // hN
                {8.2516, 2.7582},  // INCL
                {12.4768, -0.6278}  // G4
            },
            // SumStd
            {
                {2.1374, 1.6616},  // hA
                {1.2266, 4.5166},  // hN
                {2.3016, 4.5927},  // INCL
                {2.1749, 3.9374}  // G4
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
                {-0.6432, 1.9197},  // INCL
                {-0.6985, 2.2513}  // G4
            },
            // DifStd
            {
                {0.3591, 3.3994},  // hA
                {0.6809, 2.3669},  // hN
                {0.6811, 1.3115},  // INCL
                {0.3484, 1.4438}  // G4
            },
            // SumMean
            {
                {5.0295, 5.7037},  // hA
                {9.1261, 6.7006},  // hN
                {11.1969, 2.9899},  // INCL
                {13.6643, -0.3371}  // G4
            },
            // SumStd
            {
                {2.5278, 2.2597},  // hA
                {1.7786, 4.6312},  // hN
                {3.0413, 3.9273},  // INCL
                {3.2988, 3.0931}  // G4
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
    TH1D *hist = new TH1D(Form("hist_%s",label.c_str()),Form("hist_%s",label.c_str()),12, -0.2,1);
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

    double hAEst=1;
    double AltEst=1;
    if (mode==kAbs){
    const auto& tp = getPion(pdg,target);

        const LinearParams& altSumMean = (reweightModel=="hN2018") ? tp.SumMean.hN   :
                                         (reweightModel=="Geant4")  ? tp.SumMean.G4   :
                                                                       tp.SumMean.INCL;
        const LinearParams& altSumStd  = (reweightModel=="hN2018") ? tp.SumStd.hN    :
                                         (reweightModel=="Geant4")  ? tp.SumStd.G4    :
                                                                       tp.SumStd.INCL;

        double truncNormAlt = 1.0 - TMath::Freq((2.0 - meanAltSum) / stdAltSum);
        double truncNormhA  = 1.0 - TMath::Freq((2.0 - meanhASum)  / stdhASum);
        
        double hAEstSum  = TMath::Gaus(mult, meanhASum,  stdhASum,  true) / truncNormhA;
        double AltEstSum = TMath::Gaus(mult, meanAltSum, stdAltSum, true) / truncNormAlt;
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

        double truncNormAlt = 1.0 - TMath::Freq((2.0 - meanAltSum) / stdAltSum);
        double truncNormhA  = 1.0 - TMath::Freq((2.0 - meanhASum)  / stdhASum);
        
        double hAEstSum  = TMath::Gaus(mult, meanhASum,  stdhASum,  true) / truncNormhA;
        double AltEstSum = TMath::Gaus(mult, meanAltSum, stdAltSum, true) / truncNormAlt;
        double hAEstDif  = TMath::Gaus(diff, meanhADif,  stdhADif,  1);
        double AltEstDif = TMath::Gaus(diff, meanAltDif, stdAltDif, 1);
        hAEst=hAEstSum*hAEstDif; AltEst=AltEstSum*AltEstDif;
    
    }
    else{


        const auto& tp = getNucleon(pdg,target);

        const GammaParams& altGamma = (reweightModel=="hN2018") ? tp.SumGamma.hN   :
                                      (reweightModel=="Geant4")  ? tp.SumGamma.G4   :
                                                                   tp.SumGamma.INCL;

        double gammaAlt = altGamma.intercept * TMath::Exp(TMath::Power(KEini, altGamma.exp) * altGamma.slope) + altGamma.floor;
        double gammaHA  = tp.SumGamma.hA.intercept * TMath::Exp(TMath::Power(KEini, tp.SumGamma.hA.exp) * tp.SumGamma.hA.slope) + tp.SumGamma.hA.floor;
        hAEst  = TMath::Exp(-(mult - 3) * gammaHA);
        AltEst = TMath::Exp(-(mult - 3) * gammaAlt);
    }
    double rw     = (AltEst / hAEst);

    if(rw<0.001) rw=0.001; if (rw>100) rw=100;
           
    
        hist->Fill(Ebias, rw);
        } else {
            hist->Fill(Ebias, 1);
        }
    }

    return hist;
}

// =======================================================
// Main plotting function
// =======================================================
void plotVisibleEwAllRW(int probe=211, double ke=0.175, int target=1000180400,
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
    std::string particle;
    if (probe == 2212) particle = "protonPlus";
    else if (probe == 2112) particle = "neutron";
    else if (probe == 211) particle = "piPlus";
    else if (probe == 111) particle = "pi0";
    else if (probe == -211) particle = "piMinus";

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
    TH1D* h1 = makeHistogram(t1,probe,"20i",mode,useThreshold,A,target,compound,1,particle,reweightModel);
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


    h1->GetXaxis()->SetTitle("(E^{vis}-E^{ini})/E^{ini}");
    if (probe>211) h1->GetXaxis()->SetTitle("(E^{vis}-T^{ini})/T^{ini}");

    h1->GetYaxis()->SetTitle("Fraction of Absorption Int.");
    if (mode!=kAbs)    h1->GetYaxis()->SetTitle("Fraction of Knockout Int.");
    

    h1->GetXaxis()->CenterTitle();
    h1->GetYaxis()->CenterTitle();
    if (A==12) h1->GetYaxis()->SetRangeUser(0,1);
    else h1->GetYaxis()->SetRangeUser(0,0.6);
    if (h4->GetMaximum()>0.5) h1->GetYaxis()->SetRangeUser(0,1);
    h1->GetYaxis()->SetTitleOffset(1.3);

    h1->Draw("HIST ");
    h2->Draw("HIST SAME ");
    h3->Draw("HIST SAME ");
    h4->Draw("HIST SAME ");

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
leg->SetTextSize(0.035);     // smaller, fits entries


leg->AddEntry(h1,Form("hA2018 Reweighted to Mult. of %s",reweightModel.c_str()),"l");
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
    std::string channelStr="visESumThreshPionAbs";

    if (mode == kAbs)
        channelStr = useThreshold ?
            "visESumThreshPionAb" :
            "visESumPionAb";
    if (mode==kKO)
        channelStr = useThreshold ?
            "visESumThreshKO" :
            "visESumKO";
    // Argon Sum parameters




    c1->Print(Form("../plotting/%s_%s_%s_%d_%d_%dRWNuclRW%s.pdf",
                  channelStr.c_str(),
                  probeStr.c_str(),
                  keStr.c_str(),
                  target,useThreshold,compound,reweightModel.c_str()));



}