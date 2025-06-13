#ifndef POSITIONS_H
#define POSITIONS_H

#include "simparam.h"
#include "../Eigen/Sparse"

extern double TLfa[NMESHP], TLfb[NMESHP], TLdfa[NMESHP], TLdfb[NMESHP];
extern double TRfb[NMESHP], TRfc[NMESHP], TRdfb[NMESHP], TRdfc[NMESHP];
extern double DLfaba[NMESHP], DLfabb[NMESHP];
extern double DLfbba[NMESHP], DLfbbb[NMESHP];
extern double DLdfaba[NMESHP], DLdfabb[NMESHP];
extern double DLdfbba[NMESHP], DLdfbbb[NMESHP];
extern double DRfbbb[NMESHP], DRfbbc[NMESHP];
extern double DRfcbb[NMESHP], DRfcbc[NMESHP];
extern double DRdfbbb[NMESHP], DRdfbbc[NMESHP];
extern double DRdfcbb[NMESHP], DRdfcbc[NMESHP];
extern double VLfava[NMESHP], VLfavb[NMESHP];
extern double VLfbva[NMESHP], VLfbvb[NMESHP];
extern double VLdfava[NMESHP], VLdfavb[NMESHP];
extern double VLdfbva[NMESHP], VLdfbvb[NMESHP];
extern double VRfbvb[NMESHP], VRfbvc[NMESHP];
extern double VRfcvb[NMESHP], VRfcvc[NMESHP];
extern double VRdfbvb[NMESHP], VRdfbvc[NMESHP];
extern double VRdfcvb[NMESHP], VRdfcvc[NMESHP];
extern double SLsa[NMESHP], SLsb[NMESHP];
extern double SRsb[NMESHP], SRsc[NMESHP];
extern double ELfa[NMESHP], ELfb[NMESHP], ELdfa[NMESHP], ELdfb[NMESHP];
extern double ERfb[NMESHP], ERfc[NMESHP], ERdfb[NMESHP], ERdfc[NMESHP];

extern double dTLfa[NMESHP], dTLfb[NMESHP], dTLdfa[NMESHP], dTLdfb[NMESHP];
extern double dTRfb[NMESHP], dTRfc[NMESHP], dTRdfb[NMESHP], dTRdfc[NMESHP];
extern double dDLfaba[NMESHP], dDLfabb[NMESHP];
extern double dDLfbba[NMESHP], dDLfbbb[NMESHP];
extern double dDLdfaba[NMESHP], dDLdfabb[NMESHP];
extern double dDLdfbba[NMESHP], dDLdfbbb[NMESHP];
extern double dDRfbbb[NMESHP], dDRfbbc[NMESHP];
extern double dDRfcbb[NMESHP], dDRfcbc[NMESHP];
extern double dDRdfbbb[NMESHP], dDRdfbbc[NMESHP];
extern double dDRdfcbb[NMESHP], dDRdfcbc[NMESHP];
extern double dVLfava[NMESHP], dVLfavb[NMESHP];
extern double dVLfbva[NMESHP], dVLfbvb[NMESHP];
extern double dVLdfava[NMESHP], dVLdfavb[NMESHP];
extern double dVLdfbva[NMESHP], dVLdfbvb[NMESHP];
extern double dVRfbvb[NMESHP], dVRfbvc[NMESHP];
extern double dVRfcvb[NMESHP], dVRfcvc[NMESHP];
extern double dVRdfbvb[NMESHP], dVRdfbvc[NMESHP];
extern double dVRdfcvb[NMESHP], dVRdfcvc[NMESHP];
extern double dSLsa[NMESHP], dSLsb[NMESHP];
extern double dSRsb[NMESHP], dSRsc[NMESHP];
extern double dELfa[NMESHP], dELfb[NMESHP], dELdfa[NMESHP], dELdfb[NMESHP];
extern double dERfb[NMESHP], dERfc[NMESHP], dERdfb[NMESHP], dERdfc[NMESHP];

extern double Br[NMESHP];
extern double aR[NMESHP];

extern Eigen::SparseMatrix<double, Eigen::RowMajor> Ts;
extern Eigen::SparseMatrix<double, Eigen::RowMajor> Ss;
extern Eigen::SparseMatrix<double, Eigen::RowMajor> Es;

#endif // POSITIONS_H