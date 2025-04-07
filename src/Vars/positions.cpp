#include "positions.h"

double TLfa[NMESHP], TLfb[NMESHP], TLdfa[NMESHP], TLdfb[NMESHP];
double TRfb[NMESHP], TRfc[NMESHP], TRdfb[NMESHP], TRdfc[NMESHP];
double DLfaba[NMESHP], DLfabb[NMESHP];
double DLfbba[NMESHP], DLfbbb[NMESHP];
double DLdfaba[NMESHP], DLdfabb[NMESHP];
double DLdfbba[NMESHP], DLdfbbb[NMESHP];
double DRfbbb[NMESHP], DRfbbc[NMESHP];
double DRfcbb[NMESHP], DRfcbc[NMESHP];
double DRdfbbb[NMESHP], DRdfbbc[NMESHP];
double DRdfcbb[NMESHP], DRdfcbc[NMESHP];
double VLfava[NMESHP], VLfavb[NMESHP];
double VLfbva[NMESHP], VLfbvb[NMESHP];
double VLdfava[NMESHP], VLdfavb[NMESHP];
double VLdfbva[NMESHP], VLdfbvb[NMESHP];
double VRfbvb[NMESHP], VRfbvc[NMESHP];
double VRfcvb[NMESHP], VRfcvc[NMESHP];
double VRdfbvb[NMESHP], VRdfbvc[NMESHP];
double VRdfcvb[NMESHP], VRdfcvc[NMESHP];
double SLsa[NMESHP], SLsb[NMESHP];
double SRsb[NMESHP], SRsc[NMESHP];
double ELfa[NMESHP], ELfb[NMESHP], ELdfa[NMESHP], ELdfb[NMESHP];
double ERfb[NMESHP], ERfc[NMESHP], ERdfb[NMESHP], ERdfc[NMESHP];

double dTLfa[NMESHP], dTLfb[NMESHP], dTLdfa[NMESHP], dTLdfb[NMESHP];
double dTRfb[NMESHP], dTRfc[NMESHP], dTRdfb[NMESHP], dTRdfc[NMESHP];
double dDLfaba[NMESHP], dDLfabb[NMESHP];
double dDLfbba[NMESHP], dDLfbbb[NMESHP];
double dDLdfaba[NMESHP], dDLdfabb[NMESHP];
double dDLdfbba[NMESHP], dDLdfbbb[NMESHP];
double dDRfbbb[NMESHP], dDRfbbc[NMESHP];
double dDRfcbb[NMESHP], dDRfcbc[NMESHP];
double dDRdfbbb[NMESHP], dDRdfbbc[NMESHP];
double dDRdfcbb[NMESHP], dDRdfcbc[NMESHP];
double dVLfava[NMESHP], dVLfavb[NMESHP];
double dVLfbva[NMESHP], dVLfbvb[NMESHP];
double dVLdfava[NMESHP], dVLdfavb[NMESHP];
double dVLdfbva[NMESHP], dVLdfbvb[NMESHP];
double dVRfbvb[NMESHP], dVRfbvc[NMESHP];
double dVRfcvb[NMESHP], dVRfcvc[NMESHP];
double dVRdfbvb[NMESHP], dVRdfbvc[NMESHP];
double dVRdfcvb[NMESHP], dVRdfcvc[NMESHP];
double dSLsa[NMESHP], dSLsb[NMESHP];
double dSRsb[NMESHP], dSRsc[NMESHP];
double dELfa[NMESHP], dELfb[NMESHP], dELdfa[NMESHP], dELdfb[NMESHP];
double dERfb[NMESHP], dERfc[NMESHP], dERdfb[NMESHP], dERdfc[NMESHP];

double Br[NMESHP];
double aR[NMESHP];

Eigen::SparseMatrix<double, Eigen::RowMajor> Ts(2*NMESHP,2*NMESHP);
Eigen::SparseMatrix<double, Eigen::RowMajor> Ss(2*NMESHP,2*NMESHP);
Eigen::SparseMatrix<double, Eigen::RowMajor> Es(2*NMESHP,2*NMESHP);
