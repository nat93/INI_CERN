// C++
#include "iostream"
#include "string"
#include "fstream"
#include "vector"
#include "stdlib.h"
#include "time.h"
#include "iomanip"
#include "assert.h"

// ROOT
#include "TROOT.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TF1.h"
#include "TSystem.h"
#include "TLine.h"
#include "TEllipse.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TPad.h"
#include "TF1.h"
#include "TF2.h"
#include "TGraph2DErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TPolyLine.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLatex.h"
#include "TPaveStats.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TMath.h"
#include "THStack.h"
#include "TSystem.h"
#include "TBenchmark.h"
#include "TRandom3.h"
#include "TLeaf.h"
#include "TChain.h"
#include "TMultiGraph.h"
#include "TSpline.h"


// ROOT stored Variables ****************************************************************/
    //
    // Variable names are
    //      evt             -> event information
    //                              -> evt.run       run number
    //                              -> evt.evtnum    event number
    //				    -> evt.nuclear   lateral trigger event
    //      date	    -> char array of event date in DDMMYY format
    //
    //      time	    -> char array of event time in HHMMSS format
    //
    //      gPos            -> goniometer coordinates
    //                              ->  gPos.x       horizontal rotation [rad]
    //                              ->  gPos.y       tilt [rad]
    //                              ->  gPos.z       lateral translation [mm]
    //      hits            -> hit multiplicity information for planes and incoming track parameters for multi-hit (downstream) events (isHit)
    //                              ->  hits.p_nHits[plane]    number of hits on planes 0->4 (0 is upstream)
    //                              ->  hits.thetaIn_x         incoming theta in x [rad]
    //                              ->  hits.thetaIn_y         incoming theta in y [rad]
    //                              ->  hits.thetaInErr_x      incoming theta error in x [rad]
    //                              ->  hits.thetaInErr_y      incoming theta error in y [rad]
    //                              ->  hits.d0_x              impact point on xtal in x [mm]
    //                              ->  hits.d0_y              impact point on xtal in y [mm]
    //                              ->  hits.d0Err_x           impact point error on xtal in x [mm]
    //                              ->  hits.d0Err_y           impact point error on xtal in y [mm]
    //      tracks	    -> track parameters for single track events (isTrack)
    //                              ->  tracks.thetaIn_x        incoming theta in x [rad]
    //                              ->  tracks.thetaIn_y        incoming theta in y [rad]
    //                              ->  tracks.thetaOut_x       outgoing theta in x [rad]
    //                              ->  tracks.thetaOut_y       outgoing theta in y [rad]
    //                              ->  tracks.thetaInErr_x     incoming theta error in x (nominal) [rad]
    //                              ->  tracks.thetaInErr_y     incoming theta error in y (nominal) [rad]
    //                              ->  tracks.thetaOutErr_x    outgoing theta error in x (nominal) [rad]
    //                              ->  tracks.thetaOutErr_y    outgoing theta error in y (nominal) [rad]
    //                              ->  tracks.d0_x             impact point on xtal in x [mm]
    //                              ->  tracks.d0_y             impact point on xtal in y [mm]
    //                              ->  tracks.d0Err_x          impact point error on xtal in x [mm]
    //                              ->  tracks.d0Err_y          impact point error on xtal in y [mm]
    //                              ->  tracks.chi2_x           chi2 of fit in x
    //                              ->  tracks.chi2_y           chi2 of fit in y
//****************************************************************************************/

