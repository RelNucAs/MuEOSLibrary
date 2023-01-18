#pragma once
#include "parameters.hpp"
using namespace parameters;
//Abscissas and weights for Gauss-Legendre integration from x1 = 0.000 to x2 = 1.000
const double xgle[ngle] = {3.43570040745253760693880576110890438E-03,
  1.80140363610431043661669344013613777E-02,
  4.38827858743370470661237793983509225E-02,
  8.04415140888905883027354691492396756E-02,
  0.126834046769924603692847464728426333,
  0.181973159636742487273581651885776137,
  0.244566499024586450997817974522365556,
  0.313146955642290219663725911487536393,
  0.386107074429177460959751902315712705,
  0.461736739433251333122679795300580886,
  0.538263260566748666877320204699419066,
  0.613892925570822539040248097684287343,
  0.686853044357709780336274088512463607,
  0.755433500975413549002182025477634492,
  0.818026840363257512726418348114223863,
  0.873165953230075396307152535271573667,
  0.919558485911109411697264530850760276,
  0.956117214125662952933876220601649126,
  0.981985963638956895633833065598638622,
  0.996564299592547462393061194238891096};

const double wgle[ngle] = {  8.80700356957529658074597701502470771E-03,
   2.03007149001934706622614887378108318E-02,
   3.13360241670545317847510323919922958E-02,
   4.16383707883523743623790630443325598E-02,
   5.09650599086165941750439478926843639E-02,
   5.90972659807588494369806251062225373E-02,
   6.58443192245882841909956353137605861E-02,
   7.10480546591910242540764854959608998E-02,
   7.45864932363018733747738303166438139E-02,
   7.63766935653629253490395339725770993E-02,
   7.63766935653629253490395339725770993E-02,
   7.45864932363018733747738303166438139E-02,
   7.10480546591910242540764854959608998E-02,
   6.58443192245882841909956353137605861E-02,
   5.90972659807588494369806251062225373E-02,
   5.09650599086165941750439478926843639E-02,
   4.16383707883523743623790630443325598E-02,
   3.13360241670545317847510323919922958E-02,
   2.03007149001934706622614887378108318E-02,
   8.80700356957529658074597701502470771E-03};


