#include "functions.h"

void init_positions() {
    for (int i = 0; i < NMESHP; ++i) {
        // if equal position intervals
        aR[i] = -1 + (double)(i) * 2.0 / (double)(NMESHP - 1); // Added -1 to make it [-1, 1] because the loop goes from 0 to NMESHP - 1
    }
    // cout << aR[0] << " " << aR[NMESHP - 1] << endl;
    // adjust positions and intervals to get exactly [-a a] as total interval
    double scale = -a / aR[0];
    for (int i = 0; i < NMESHP; ++i) { //	aR  = (-a+2*a/N/2:2*a/N:a-2*a/N/2)'; % cm
        aR[i] = aR[i] * scale;
    }
    for (int i = 0; i < NMESHP; ++i) { //	% cm
        aR[i] = aR[i] + R;
    }

    // double Br[NMESHP];
    for (int id = 0; id < NMESHP; ++id) {
        Br[id] = Bt * R / (aR[id]);
    }

    initAuxiliarPos();

    for (int i = 1; i < NMESHP - 1; ++i) { // row
        for (int j = 0; j < NMESHP; ++j) { // column
            if (i == j) {
                Ts.insert(2 * i, 2 * j) = TLfb[i] + TRfb[i];
                Ts.insert(2 * i + 1, 2 * j) = dTLfb[i] + dTRfb[i];
                Ts.insert(2 * i, 2 * j + 1) = TLdfb[i] + TRdfb[i];
                Ts.insert(2 * i + 1, 2 * j + 1) = dTLdfb[i] + dTRdfb[i];
            } else if (j == i - 1) {
                Ts.insert(2 * i, 2 * j) = TLfa[i];
                Ts.insert(2 * i + 1, 2 * j) = dTLfa[i];
                Ts.insert(2 * i, 2 * j + 1) = TLdfa[i];
                Ts.insert(2 * i + 1, 2 * j + 1) = dTLdfa[i];
            } else if (j == i + 1) {
                Ts.insert(2 * i, 2 * j) = TRfc[i];
                Ts.insert(2 * i + 1, 2 * j) = dTRfc[i];
                Ts.insert(2 * i, 2 * j + 1) = TRdfc[i];
                Ts.insert(2 * i + 1, 2 * j + 1) = dTRdfc[i];
            }
        }
    }
    for (int i = 0; i < 1; ++i) {     // row
        for (int j = 0; j < 2; ++j) { // column
            if (i == j) {
                Ts.insert(2 * i, 2 * j) = TRfb[i];
                Ts.insert(2 * i + 1, 2 * j) = dTRfb[i];
                Ts.insert(2 * i, 2 * j + 1) = TRdfb[i];
                Ts.insert(2 * i + 1, 2 * j + 1) = dTRdfb[i];
            } else if (j == i + 1) {
                Ts.insert(2 * i, 2 * j) = TRfc[i];
                Ts.insert(2 * i + 1, 2 * j) = dTRfc[i];
                Ts.insert(2 * i, 2 * j + 1) = TRdfc[i];
                Ts.insert(2 * i + 1, 2 * j + 1) = dTRdfc[i];
            }
        }
    }
    for (int i = NMESHP - 1; i < NMESHP; ++i) {     // row
        for (int j = NMESHP - 2; j < NMESHP; ++j) { // column
            if (i == j) {
                Ts.insert(2 * i, 2 * j) = TLfb[i];
                Ts.insert(2 * i + 1, 2 * j) = dTLfb[i];
                Ts.insert(2 * i, 2 * j + 1) = TLdfb[i];
                Ts.insert(2 * i + 1, 2 * j + 1) = dTLdfb[i];
            } else if (j == i - 1) {
                Ts.insert(2 * i, 2 * j) = TLfa[i];
                Ts.insert(2 * i + 1, 2 * j) = dTLfa[i];
                Ts.insert(2 * i, 2 * j + 1) = TLdfa[i];
                Ts.insert(2 * i + 1, 2 * j + 1) = dTLdfa[i];
            }
        }
    }

    for (int i = 1; i < NMESHP - 1; ++i) { // row
        for (int j = 0; j < NMESHP; ++j) { // column
            if (i == j) {
                Ss.insert(2 * i, 2 * j) = SLsb[i] + SRsb[i];
                Ss.insert(2 * i + 1, 2 * j) = dSLsb[i] + dSRsb[i];
            } else if (j == i - 1) {
                Ss.insert(2 * i, 2 * j) = SLsa[i];
                Ss.insert(2 * i + 1, 2 * j) = dSLsa[i];
            } else if (j == i + 1) {
                Ss.insert(2 * i, 2 * j) = SRsc[i];
                Ss.insert(2 * i + 1, 2 * j) = dSRsc[i];
            }
        }
    }
    for (int i = 0; i < 1; ++i) {     // row
        for (int j = 0; j < 2; ++j) { // column
            if (i == j) {
                Ss.insert(2 * i, 2 * j) = SRsb[i];
                Ss.insert(2 * i + 1, 2 * j) = dSRsb[i];
            } else if (j == i + 1) {
                Ss.insert(2 * i, 2 * j) = SRsc[i];
                Ss.insert(2 * i + 1, 2 * j) = dSRsc[i];
            }
        }
    }
    for (int i = NMESHP - 1; i < NMESHP; ++i) {     // row
        for (int j = NMESHP - 2; j < NMESHP; ++j) { // column
            if (i == j) {
                Ss.insert(2 * i, 2 * j) = SLsb[i];
                Ss.insert(2 * i + 1, 2 * j) = dSLsb[i];
            } else if (j == i - 1) {
                Ss.insert(2 * i, 2 * j) = SLsa[i];
                Ss.insert(2 * i + 1, 2 * j) = dSLsa[i];
            }
        }
    }
}

void initAuxiliarPos() {
    double r_a, r_b, r_c;

    double t1, t2; //, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t30, t31, t32, t33, t34, t35, t36, t37, t38, t39, t40, t41, t42, t43, t44, t45, t46, t47, t48, t49, t50, t51, t52, t53, t54, t55, t56, t57, t58, t59, t60, t61, t62, t63, t64, t65, t66, t67, t68, t69, t70, t71, t72, t73, t74, t75, t76, t77, t78, t79, t80, t81, t82, t83, t84, t85, t86, t87, t88, t89, t90, t91, t92, t93, t94, t95, t96, t97, t98, t99, t100, t101, t102, t103, t104, t105, t106, t107, t108, t109, t110, t111, t112, t113, t114, t115, t116, t117, t118, t119, t120, t121, t122, t123, t124, t125, t126, t127, t128, t129, t130, t131, t132, t133, t134, t135, t136, t137, t138, t139, t140, t141, t142, t143, t144, t145, t146, t147, t148, t149, t150;

    // left sides
    for (int i = 1; i < NMESHP; ++i) {

        r_a = aR[i - 1];
        r_b = aR[i];

        t1 = -0.9e1 / 0.140e3;
        TLfa[i] = t1 * (r_a * r_a - r_b * r_b);
        t1 = 0.1e1 / 0.35e2;
        TLfb[i] = t1 * (r_b - r_a) * (0.3e1 * r_a + 0.10e2 * r_b);
        t1 = -0.1e1 / 0.420e3;
        TLdfa[i] = t1 * (r_b - r_a) * (0.6e1 * r_a * r_a + (r_a - 0.7e1 * r_b) * r_b);
        t1 = 0.1e1 / 0.420e3;
        TLdfb[i] = t1 * (r_b - r_a) * (0.7e1 * r_a * r_a + (0.8e1 * r_a - 0.15e2 * r_b) * r_b);
        t1 = r_b - r_a;
        t1 = 0.1e1 / t1;
        DLfaba[i] = (0.12e2 / 0.35e2 * r_a + 0.9e1 / 0.35e2 * r_b) * t1;
        t1 = r_b - r_a;
        t1 = 0.1e1 / t1;
        DLfabb[i] = (0.9e1 / 0.35e2 * r_a + 0.12e2 / 0.35e2 * r_b) * t1;
        t1 = r_b - r_a;
        t1 = 0.1e1 / t1;
        DLfbba[i] = (-0.12e2 / 0.35e2 * r_a - 0.9e1 / 0.35e2 * r_b) * t1;
        t1 = r_b - r_a;
        t1 = 0.1e1 / t1;
        DLfbbb[i] = (-0.9e1 / 0.35e2 * r_a - 0.12e2 / 0.35e2 * r_b) * t1;
        t1 = 0.1e1 / 0.35e2;
        DLdfaba[i] = t1 * (r_b - r_a);
        DLdfabb[i] = r_a / 0.35e2 + r_b / 0.14e2;
        DLdfbba[i] = r_a / 0.14e2 + r_b / 0.35e2;
        DLdfbbb[i] = 0.34e2 / 0.35e2 * r_b + r_a / 0.35e2;
        VLfava[i] = -0.3e1 / 0.14e2 * r_a - r_b / 0.10e2;
        VLfavb[i] = -r_a / 0.10e2 - 0.3e1 / 0.35e2 * r_b;
        VLfbva[i] = -0.3e1 / 0.35e2 * r_a - r_b / 0.10e2;
        VLfbvb[i] = 0.11e2 / 0.14e2 * r_b - r_a / 0.10e2;
        t1 = -0.1e1 / 0.140e3;
        VLdfava[i] = t1 * (r_b - r_a) * (0.5e1 * r_a + 0.3e1 * r_b);
        t1 = 0.3e1 / 0.140e3;
        VLdfavb[i] = t1 * (r_a * r_a - r_b * r_b);
        t1 = -0.3e1 / 0.140e3;
        VLdfbva[i] = t1 * (r_a * r_a - r_b * r_b);
        t1 = 0.1e1 / 0.140e3;
        VLdfbvb[i] = t1 * (r_b - r_a) * (0.3e1 * r_a + 0.5e1 * r_b);
        SLsa[i] = (r_a / 0.15e2 + r_b / 0.12e2) * (r_b - r_a);
        SLsb[i] = (r_a / 0.12e2 + 0.4e1 / 0.15e2 * r_b) * (r_b - r_a);
        //ELfa[i] = -0.11e2 / 0.35e2 * r_a - 0.13e2 / 0.70e2 * r_b; //assumption: not used
        //ELfb[i] = 0.24e2 / 0.35e2 * r_b - 0.13e2 / 0.70e2 * r_a;
        //ELdfa[i] = (-r_a / 0.70e2 - 0.3e1 / 0.70e2 * r_b) * r_b + 0.2e1 / 0.35e2 * r_a * r_a;
        //ELdfb[i] = (0.2e1 / 0.35e2 * r_b - r_a / 0.70e2) * r_b - 0.3e1 / 0.70e2 * r_a * r_a;
        t1 = r_b - r_a;
        dTLfa[i] = -pow(t1, 0.2e1) * (0.7e1 * r_a + 0.6e1 * r_b) / 0.420e3;
        t1 = r_b - r_a;
        dTLfb[i] = -pow(t1, 0.2e1) * (0.7e1 * r_a + 0.15e2 * r_b) / 0.420e3;
        t1 = r_b - r_a;
        dTLdfa[i] = pow(t1, 0.2e1) * (r_a * r_a - r_b * r_b) / 0.280e3;
        t1 = r_b - r_a;
        dTLdfb[i] = -pow(t1, 0.2e1) * (0.3e1 * r_a * r_a + (0.2e1 * r_a - 0.5e1 * r_b) * r_b) / 0.840e3;
        dDLfaba[i] = -r_a / 0.14e2 - r_b / 0.35e2;
        t1 = 0.1e1 / 0.35e2;
        dDLfabb[i] = t1 * (r_b - r_a);
        dDLfbba[i] = r_a / 0.14e2 + r_b / 0.35e2;
        t1 = -0.1e1 / 0.35e2;
        dDLfbbb[i] = t1 * (r_b - r_a);
        dDLdfaba[i] = (r_b / 0.420e3 + r_a / 0.84e2) * r_b - r_a * r_a / 0.70e2;
        dDLdfabb[i] = (r_b / 0.70e2 - r_a / 0.84e2) * r_b - r_a * r_a / 0.420e3;
        dDLdfbba[i] = (-r_a / 0.210e3 - r_b / 0.70e2) * r_b + 0.2e1 / 0.105e3 * r_a * r_a;
        dDLdfbbb[i] = (r_a / 0.14e2 - 0.3e1 / 0.35e2 * r_b) * r_b + r_a * r_a / 0.70e2;
        t1 = 0.1e1 / 0.420e3;
        dVLfava[i] = t1 * (r_b - r_a) * (0.23e2 * r_a + 0.8e1 * r_b);
        t1 = 0.1e1 / 0.420e3;
        dVLfavb[i] = t1 * (r_b - r_a) * (0.8e1 * r_a + 0.3e1 * r_b);
        t1 = 0.5e1;
        t2 = 0.1e1 / 0.420e3;
        dVLfbva[i] = t2 * (r_b - r_a) * (r_a * t1 - r_b);
        t1 = -0.1e1 / 0.420e3;
        dVLfbvb[i] = t1 * (r_b - r_a) * (r_a + 0.45e2 * r_b);
        t1 = r_b - r_a;
        dVLdfava[i] = pow(t1, 0.2e1) * (0.7e1 * r_a + 0.3e1 * r_b) / 0.840e3;
        t1 = r_b - r_a;
        dVLdfavb[i] = pow(t1, 0.2e1) * (0.3e1 * r_a + r_b) / 0.840e3;
        t1 = r_b - r_a;
        dVLdfbva[i] = -pow(t1, 0.2e1) * (0.3e1 * r_a + r_b) / 0.840e3;
        t1 = r_b - r_a;
        dVLdfbvb[i] = -pow(t1, 0.2e1) * (r_a - 0.5e1 * r_b) / 0.840e3;
        t1 = r_b - r_a;
        t2 = -0.1e1 / 0.60e2;
        dSLsa[i] = t2 * (r_a + r_b) * pow(t1, 0.2e1);
        t1 = r_b - r_a;
        dSLsb[i] = (-r_a / 0.60e2 - r_b / 0.30e2) * pow(t1, 0.2e1);
    }

    // right sides
    for (int i = 0; i < NMESHP - 1; ++i) {

        r_b = aR[i];
        r_c = aR[i + 1];

        t1 = 0.1e1 / 0.35e2;
        TRfb[i] = t1 * (r_c - r_b) * (0.10e2 * r_b + 0.3e1 * r_c);
        t1 = 0.9e1 / 0.140e3;
        TRfc[i] = t1 * (-r_b * r_b + r_c * r_c);
        t1 = -0.1e1 / 0.420e3;
        TRdfb[i] = t1 * (r_c - r_b) * (0.15e2 * r_b * r_b + (-0.8e1 * r_b - 0.7e1 * r_c) * r_c);
        t1 = -0.6e1;
        t2 = 0.1e1 / 0.420e3;
        TRdfc[i] = t2 * (r_c - r_b) * (0.7e1 * r_b * r_b + (r_c * t1 - r_b) * r_c);
        t1 = r_c - r_b;
        t1 = 0.1e1 / t1;
        DRfbbb[i] = (-0.12e2 / 0.35e2 * r_b - 0.9e1 / 0.35e2 * r_c) * t1;
        t1 = r_c - r_b;
        t1 = 0.1e1 / t1;
        DRfbbc[i] = (-0.9e1 / 0.35e2 * r_b - 0.12e2 / 0.35e2 * r_c) * t1;
        t1 = r_c - r_b;
        t1 = 0.1e1 / t1;
        DRfcbb[i] = (0.12e2 / 0.35e2 * r_b + 0.9e1 / 0.35e2 * r_c) * t1;
        t1 = r_c - r_b;
        t1 = 0.1e1 / t1;
        DRfcbc[i] = (0.9e1 / 0.35e2 * r_b + 0.12e2 / 0.35e2 * r_c) * t1;
        DRdfbbb[i] = -0.34e2 / 0.35e2 * r_b - r_c / 0.35e2;
        DRdfbbc[i] = -r_b / 0.35e2 - r_c / 0.14e2;
        DRdfcbb[i] = -r_b / 0.14e2 - r_c / 0.35e2;
        t1 = 0.1e1 / 0.35e2;
        DRdfcbc[i] = t1 * (r_c - r_b);
        VRfbvb[i] = -0.11e2 / 0.14e2 * r_b + r_c / 0.10e2;
        VRfbvc[i] = r_b / 0.10e2 + 0.3e1 / 0.35e2 * r_c;
        VRfcvb[i] = 0.3e1 / 0.35e2 * r_b + r_c / 0.10e2;
        VRfcvc[i] = r_b / 0.10e2 + 0.3e1 / 0.14e2 * r_c;
        VRdfbvb[i] = (0.3e1 / 0.140e3 * r_c + r_b / 0.70e2) * r_c - r_b * r_b / 0.28e2;
        t1 = 0.3e1 / 0.140e3;
        VRdfbvc[i] = t1 * (-r_b * r_b + r_c * r_c);
        t1 = -0.3e1 / 0.140e3;
        VRdfcvb[i] = t1 * (-r_b * r_b + r_c * r_c);
        VRdfcvc[i] = (r_b / 0.70e2 - r_c / 0.28e2) * r_c + 0.3e1 / 0.140e3 * r_b * r_b;
        SRsb[i] = (0.4e1 / 0.15e2 * r_b + r_c / 0.12e2) * (r_c - r_b);
        SRsc[i] = (r_b / 0.12e2 + r_c / 0.15e2) * (r_c - r_b);
        t1 = r_c - r_b;
        dTRfb[i] = (r_b / 0.28e2 + r_c / 0.60e2) * pow(t1, 0.2e1);
        dTRfc[i] = (r_b / 0.70e2 + r_c / 0.60e2) * pow(t1, 0.2e1);
        dTRdfb[i] = (-r_b * r_b / 0.168e3 - (-0.2e1 * r_b - 0.3e1 * r_c) * r_c / 0.840e3) * pow(t1, 0.2e1);
        t2 = -0.1e1 / 0.280e3;
        dTRdfc[i] = t2 * (-r_b * r_b + r_c * r_c) * pow(t1, 0.2e1);
        t1 = -0.1e1 / 0.35e2;
        dDRfbbb[i] = t1 * (r_c - r_b);
        dDRfbbc[i] = -r_b / 0.35e2 - r_c / 0.14e2;
        t1 = 0.1e1 / 0.35e2;
        dDRfcbb[i] = t1 * (r_c - r_b);
        dDRfcbc[i] = r_b / 0.35e2 + r_c / 0.14e2;
        dDRdfbbb[i] = (-r_c / 0.70e2 - r_b / 0.14e2) * r_c + 0.3e1 / 0.35e2 * r_b * r_b;
        dDRdfbbc[i] = (r_b / 0.210e3 - 0.2e1 / 0.105e3 * r_c) * r_c + r_b * r_b / 0.70e2;
        dDRdfcbb[i] = (r_c / 0.420e3 + r_b / 0.84e2) * r_c - r_b * r_b / 0.70e2;
        dDRdfcbc[i] = (r_c / 0.70e2 - r_b / 0.84e2) * r_c - r_b * r_b / 0.420e3;
        dVRfbvb[i] = (-r_c / 0.420e3 - 0.11e2 / 0.105e3 * r_b) * r_c + 0.3e1 / 0.28e2 * r_b * r_b;
        dVRfbvc[i] = (r_c / 0.84e2 - r_b / 0.70e2) * r_c + r_b * r_b / 0.420e3;
        dVRfcvb[i] = (0.2e1 / 0.105e3 * r_c - r_b / 0.84e2) * r_c - r_b * r_b / 0.140e3;
        dVRfcvc[i] = (0.23e2 / 0.420e3 * r_c - r_b / 0.28e2) * r_c - 0.2e1 / 0.105e3 * r_b * r_b;
        t1 = r_b * r_b;
        dVRdfbvb[i] = ((r_c / 0.840e3 - r_b / 0.120e3) * r_c + 0.11e2 / 0.840e3 * t1) * r_c - r_b * t1 / 0.168e3;
        dVRdfbvc[i] = r_b * r_b * (r_b + r_c) / 0.840e3 + (r_c / 0.280e3 - r_b / 0.168e3) * r_c * r_c;
        dVRdfcvb[i] = -r_b * r_b * (r_b + r_c) / 0.840e3 + (r_b / 0.168e3 - r_c / 0.280e3) * r_c * r_c;
        t1 = r_b * r_b;
        dVRdfcvc[i] = ((0.11e2 / 0.840e3 * r_b - r_c / 0.120e3) * r_c - t1 / 0.840e3) * r_c - r_b * t1 / 0.280e3;
        t1 = r_c - r_b;
        dSRsb[i] = pow(t1, 0.2e1) * (r_c + 0.2e1 * r_b) / 0.60e2;
        t1 = r_c - r_b;
        dSRsc[i] = pow(t1, 0.2e1) * (r_b + r_c) / 0.60e2;
    }
}

double computeN(double p, double Ta) {
    return p * 100.0 / (Ta * 11600.0 * kb) / 1.0e6;
}

void initializeSpecies() {
    for (int im = 0; im < NMESHP; ++im) {
        double fct_n0, fct_dn0;

        if (aR[im] <= rmaxini) {
            fct_n0 = nebackgroundl + exp(-pow((aR[im] - rmaxini) / (widthini), 2.0)); //spiked at rmaxini width widthini
            fct_dn0 = -2.0 * (aR[im] - rmaxini) / pow((widthini), 2.0) * (fct_n0 - nebackgroundl);
        } else {
            fct_n0 = nebackgroundr + exp(-pow((aR[im] - rmaxini) / (2.0 * widthini), 2.0));
            fct_dn0 = -2.0 * (aR[im] - rmaxini) / pow((2.0 * widthini), 2.0) * (fct_n0 - nebackgroundr);
        }

        Er.Ee[im] = ENERGY_FACTOR * Te0 * nr.ne[im];
        nr.nH[im] = nH0 * fct_n0;
        if (nr.nH[im] < 0.0)
            cout << "nH<0" << endl;
        nr.nHi[im] = nHi0 * fct_n0;
        nr.nH2i[im] = nH2i0 * fct_n0;
        nr.nH3i[im] = nH3i0 * fct_n0;
        nr.nH2[im] = nH20; // - 0.5*(nr.nHi[im]+nr.nH[im]) - nr.nH2i[im] - 2.0/3.0*nr.nH3i[im];//nr.nH+
        nr.nHeII[im] = nHeII0 * fct_n0;
        nr.nHeIII[im] = nHeIII0 * fct_n0;
        nr.nHeI[im] = nHeI0; //  - nr.nHeIII[im] - nr.nHeII[im];
        nr.nCII[im] = nCII0 * fct_n0;
        nr.nCIII[im] = nCIII0 * fct_n0;
        nr.nCIV[im] = nCIV0 * fct_n0;
        nr.nCV[im] = nCV0 * fct_n0;
        nr.nCI[im] = nCI0 - nr.nCII[im] - nr.nCIII[im] - nr.nCIV[im] - nr.nCV[im];
        nr.ne[im] = nr.nHi[im] + nr.nH2i[im] + nr.nH3i[im] + nr.nHeII[im] + 2.0 * nr.nHeIII[im];
        if (bimpur) {
            nr.ne[im] += nr.nCII[im] + 2.0 * nr.nCIII[im] + 3.0 * nr.nCIV[im] + 4.0 * nr.nCV[im];
        }
        nr.xnH[im] = nH0 * fct_dn0;
        nr.xnHi[im] = nHi0 * fct_dn0;
        nr.xnH2i[im] = nH2i0 * fct_dn0;
        nr.xnH3i[im] = nH3i0 * fct_dn0;
        nr.xnH2[im] = 0.0; //- 0.5*(nr.xnHi[im]) - nr.xnH2i[im] - 2.0/3.0*nr.xnH3i[im];//nr.nH+
        nr.xnHeII[im] = nHeII0 * fct_dn0;
        nr.xnHeIII[im] = nHeIII0 * fct_dn0;
        nr.xnHeI[im] = 0.0; //- nr.xnHeIII[im] - nr.xnHeII[im];
        nr.xnCII[im] = nCII0 * fct_dn0;
        nr.xnCIII[im] = nCIII0 * fct_dn0;
        nr.xnCI[im] = -nr.xnCII[im] - nr.xnCIII[im] - nr.xnCIV[im] - nr.xnCV[im];
        nr.xne[im] = nr.xnHi[im] + nr.xnH2i[im] + nr.xnH3i[im] + nr.xnHeII[im] + 2.0 * nr.xnHeIII[im];
        if (bimpur) {
            nr.xne[im] += nr.xnCII[im] + 2.0 * nr.xnCIII[im] + 3.0 * nr.xnCIV[im] + 4.0 * nr.xnCV[im];
        }
        Er.Ee[im] = ENERGY_FACTOR * Te0 * nr.ne[im];
        Er.EH[im] = ENERGY_FACTOR * 2.0 * nr.nH[im];
        Er.EH2[im] = ENERGY_FACTOR * Ta0 * nr.nH2[im];
        Er.EHi[im] = ENERGY_FACTOR * 1.0 * nr.nHi[im];
        Er.EH2i[im] = ENERGY_FACTOR * 1.0 * nr.nH2i[im];
        Er.EH3i[im] = ENERGY_FACTOR * 1.0 * nr.nH3i[im];
        Er.EHeI[im] = ENERGY_FACTOR * Ta0 * nr.nHeI[im];
        Er.EHeII[im] = ENERGY_FACTOR * 1.0 * nr.nHeII[im];
        Er.EHeIII[im] = ENERGY_FACTOR * 1.0 * nr.nHeIII[im];

        Er.xEe[im] = ENERGY_FACTOR * Te0 * nr.xne[im];
        Er.xEH[im] = ENERGY_FACTOR * 3.0 * nr.xnH[im];
        Er.xEH2[im] = ENERGY_FACTOR * Ta0 * nr.xnH2[im];
        Er.xEHi[im] = ENERGY_FACTOR * 1.0 * nr.xnHi[im];
        Er.xEH2i[im] = ENERGY_FACTOR * 1.0 * nr.xnH2i[im];
        Er.xEH3i[im] = ENERGY_FACTOR * 1.0 * nr.xnH3i[im];
        Er.xEHeI[im] = ENERGY_FACTOR * Ta0 * nr.xnHeI[im];
        Er.xEHeII[im] = ENERGY_FACTOR * 1.0 * nr.xnHeII[im];
        Er.xEHeIII[im] = ENERGY_FACTOR * 1.0 * nr.xnHeIII[im];

        Power.Pabs1[im] = 0;
        Power.Pabs2[im] = 0;
        Power.Pabs3[im] = 0;
        Power.Pabs4[im] = 0;

        PRFe_id[im] = 0.0;
        PRFHi_id[im] = 0.0;
        PRFH2i_id[im] = 0.0;
        PRFH3i_id[im] = 0.0;
        PRFHeII_id[im] = 0.0;
        PRFHeIII_id[im] = 0.0;
    }
    return;
}

double ndamp(double dens) {
    return pow(1.0 + nevac / dens, -0.66);
}

void limiters() {

    double Tauion;
    double gEe = 1.0;
    double Z = 1.0;  // charge state
    double mu = 1.0; // mass unit

    int im = 0;
    bool add = 1;
    for (int ig = 0; ig < 2; ++ig) {
        while ((aR[im] < lHFS) | (aR[im] > lLFS)) {
            Z = 1.0;
            mu = 1.0;
            Tauion = 0.66 * 2.0 * pi * aR[im] / nlimiters / (9.79e5 * sqrt(Z / mu * (Tr.Te[im] + Tr.THi[im] * nr.nHi[im] / nr.ne[im])) * pow(1.0 + nevac / nr.ne[im], -0.66) * pow(1.0 + nevac / nr.nHi[im], -0.66)); //;        if (nr.ne[im]<1.0e4) {Tauion=Tauion*1.0e4/nr.ne[im];}
            // Tauion = 0.66*2.0*pi*aR[im]/nlimiters/(9.79e5*sqrt(Z/mu*(Tr.Te[im]+Tr.THi[im])) * pow(1.0+nevac/nr.ne[im],-0.66) * pow(1.0+nevac/nr.nHi[im],-0.66)) ; //;        if (nr.ne[im]<1.0e4) {Tauion=Tauion*1.0e4/nr.ne[im];}
            dnr.dnHi[im] -= nr.nHi[im] / Tauion;
	    #ifdef debug
	    	cout << "H+ losses = " << dnr.dnHi[im] << " [1/s]" << endl;
	    #endif
            dEr.dEHi[im] -= Er.EHi[im] / (Tauion / gEe);
            dnr.dne[im] -= Z * nr.nHi[im] / Tauion;
	    #ifdef debug
	    	cout << "e- losses = " << dnr.dne[im] << " [1/s]" << endl;
	    #endif
            dEr.dEe[im] -= Z * nr.nHi[im] / (Tauion / gEe) * Tr.Te[im] * 3.0 / 2.0;
            dnr.dnH2[im] += 0.5 * nr.nHi[im] / Tauion;
            dEr.dEH2[im] += 0.5 * nr.nHi[im] / Tauion * 3.0 / 2.0 * Ta0;
            // Z=1.0;
            mu = 2.0;
            Tauion = 0.66 * 2.0 * pi * aR[im] / nlimiters / (9.79e5 * sqrt(Z / mu * (Tr.Te[im] + Tr.TH2i[im] * nr.nH2i[im] / nr.ne[im])) * pow(1.0 + nevac / nr.ne[im], -0.66) * pow(1.0 + nevac / nr.nH2i[im], -0.66)); //;        if (nr.ne[im]<1.0e4) {Tauion=Tauion*1.0e4/nr.ne[im];}
            // Tauion = 0.66*2.0*pi*aR[im]/nlimiters/(9.79e5*sqrt(Z/mu*(Tr.Te[im]+Tr.TH2i[im])) * pow(1.0+nevac/nr.ne[im],-0.66) * pow(1.0+nevac/nr.nH2i[im],-0.66)) ; //;        if (nr.ne[im]<1.0e4) {Tauion=Tauion*1.0e4/nr.ne[im];}
            dnr.dnH2i[im] -= nr.nH2i[im] / Tauion;
            dEr.dEH2i[im] -= Er.EH2i[im] / (Tauion / gEe);
            dnr.dne[im] -= Z * nr.nH2i[im] / Tauion;
            dEr.dEe[im] -= Z * nr.nH2i[im] / (Tauion / gEe) * Tr.Te[im] * 3.0 / 2.0;
            dnr.dnH2[im] += 1.0 * nr.nH2i[im] / Tauion;
            dEr.dEH2[im] += 1.0 * nr.nH2i[im] / Tauion * 3.0 / 2.0 * Ta0;
            // Z=1.0;
            mu = 3.0;
            Tauion = 0.66 * 2.0 * pi * aR[im] / nlimiters / (9.79e5 * sqrt(Z / mu * (Tr.Te[im] + Tr.TH3i[im] * nr.nH3i[im] / nr.ne[im])) * pow(1.0 + nevac / nr.ne[im], -0.66) * pow(1.0 + nevac / nr.nH3i[im], -0.66)); //;        if (nr.ne[im]<1.0e4) {Tauion=Tauion*1.0e4/nr.ne[im];}
            // Tauion = 0.66*2.0*pi*aR[im]/nlimiters/(9.79e5*sqrt(Z/mu*(Tr.Te[im]+Tr.TH3i[im])) * pow(1.0+nevac/nr.ne[im],-0.66) * pow(1.0+nevac/nr.nH3i[im],-0.66)) ; //;        if (nr.ne[im]<1.0e4) {Tauion=Tauion*1.0e4/nr.ne[im];}
            dnr.dnH3i[im] -= nr.nH3i[im] / Tauion;
            dEr.dEH3i[im] -= Er.EH3i[im] / (Tauion / gEe);
            dnr.dne[im] -= Z * nr.nH3i[im] / Tauion;
            dEr.dEe[im] -= Z * nr.nH3i[im] / (Tauion / gEe) * Tr.Te[im] * 3.0 / 2.0;
            dnr.dnH2[im] += 1.5 * nr.nH3i[im] / Tauion;
            dEr.dEH2[im] += 1.5 * nr.nH3i[im] / Tauion * 3.0 / 2.0 * Ta0;
            // Z=1.0;
            mu = 4.0;
            Tauion = 0.66 * 2.0 * pi * aR[im] / nlimiters / (9.79e5 * sqrt(Z / mu * (Tr.Te[im] + Tr.THeII[im] * nr.nHeII[im] / nr.ne[im])) * pow(1.0 + nevac / nr.ne[im], -0.66) * pow(1.0 + nevac / nr.nHeII[im], -0.66)); //;        if (nr.ne[im]<1.0e4) {Tauion=Tauion*1.0e4/nr.ne[im];}
            // Tauion = 0.66*2.0*pi*aR[im]/nlimiters/(9.79e5*sqrt(Z/mu*(Tr.Te[im]+Tr.THeII[im])) * pow(1.0+nevac/nr.ne[im],-0.66) * pow(1.0+nevac/nr.nHeII[im],-0.66)) ; //;        if (nr.ne[im]<1.0e4) {Tauion=Tauion*1.0e4/nr.ne[im];}
            dnr.dnHeII[im] -= nr.nHeII[im] / Tauion;
            dEr.dEHeII[im] -= Er.EHeII[im] / (Tauion / gEe);
            dnr.dne[im] -= Z * nr.nHeII[im] / Tauion;
            dEr.dEe[im] -= Z * nr.nHeII[im] / (Tauion / gEe) * Tr.Te[im] * 3.0 / 2.0;
            dnr.dnHeI[im] += nr.nHeII[im] / Tauion;
            dEr.dEHeI[im] += nr.nHeII[im] / Tauion * 3.0 / 2.0 * Ta0;
            Z = 2.0;
            // mu=4.0;
            Tauion = 0.66 * 2.0 * pi * aR[im] / nlimiters / (9.79e5 * sqrt(Z / mu * (Tr.Te[im] + Tr.THeIII[im] * nr.nHeIII[im] / nr.ne[im])) * pow(1.0 + nevac / nr.ne[im], -0.66) * pow(1.0 + nevac / nr.nHeIII[im], -0.66)); //;        if (nr.ne[im]<1.0e4) {Tauion=Tauion*1.0e4/nr.ne[im];}
            // Tauion = 0.66*2.0*pi*aR[im]/nlimiters/(9.79e5*sqrt(Z/mu*(Tr.Te[im]+Tr.THeIII[im])) * pow(1.0+nevac/nr.ne[im],-0.66) * pow(1.0+nevac/nr.nHeIII[im],-0.66)) ; //;        if (nr.ne[im]<1.0e4) {Tauion=Tauion*1.0e4/nr.ne[im];}
            dnr.dnHeIII[im] -= nr.nHeIII[im] / Tauion;
            dEr.dEHeIII[im] -= Er.EHeIII[im] / (Tauion / gEe);
            dnr.dne[im] -= Z * nr.nHeIII[im] / Tauion;
            dEr.dEe[im] -= Z * nr.nHeIII[im] / (Tauion / gEe) * Tr.Te[im] * 3.0 / 2.0;
            dnr.dnHeI[im] += nr.nHeIII[im] / Tauion;
            dEr.dEHeI[im] += nr.nHeIII[im] / Tauion * 3.0 / 2.0 * Ta0;
            if (bimpur) {
                Z = 1.0;
                mu = 12.0;
                Tauion = 0.66 * 2.0 * pi * aR[im] / nlimiters / (9.79e5 * sqrt(Z / mu * Tr.Te[im]) * pow(1.0 + nevac / nr.ne[im], -0.66) * pow(1.0 + nevac / nr.nCII[im], -0.66)); //;        if (nr.ne[im]<1.0e4) {Tauion=Tauion*1.0e4/nr.ne[im];}
                dnr.dnCII[im] -= nr.nCII[im] / Tauion;
                dnr.dne[im] -= Z * nr.nCII[im] / Tauion;
                dEr.dEe[im] -= Z * nr.nCII[im] / (Tauion / gEe) * Tr.Te[im] * 3.0 / 2.0;
                dnr.dnCI[im] += nr.nCII[im] / Tauion;
                Z = 2.0;
                // mu=12.0;
                Tauion = 0.66 * 2.0 * pi * aR[im] / nlimiters / (9.79e5 * sqrt(Z / mu * Tr.Te[im]) * pow(1.0 + nevac / nr.ne[im], -0.66) * pow(1.0 + nevac / nr.nCIII[im], -0.66)); //;        if (nr.ne[im]<1.0e4) {Tauion=Tauion*1.0e4/nr.ne[im];}
                dnr.dnCIII[im] -= nr.nCIII[im] / Tauion;
                dnr.dne[im] -= Z * nr.nCIII[im] / Tauion;
                dEr.dEe[im] -= Z * nr.nCIII[im] / (Tauion / gEe) * Tr.Te[im] * 3.0 / 2.0;
                dnr.dnCI[im] += nr.nCIII[im] / Tauion;
                Z = 3.0;
                // mu=12.0;
                Tauion = 0.66 * 2.0 * pi * aR[im] / nlimiters / (9.79e5 * sqrt(Z / mu * Tr.Te[im]) * pow(1.0 + nevac / nr.ne[im], -0.66) * pow(1.0 + nevac / nr.nCIV[im], -0.66)); //;        if (nr.ne[im]<1.0e4) {Tauion=Tauion*1.0e4/nr.ne[im];}
                dnr.dnCIV[im] -= nr.nCIV[im] / Tauion;
                dnr.dne[im] -= Z * nr.nCIV[im] / Tauion;
                dEr.dEe[im] -= Z * nr.nCIV[im] / (Tauion / gEe) * Tr.Te[im] * 3.0 / 2.0;
                dnr.dnCI[im] += nr.nCIV[im] / Tauion;
                Z = 4.0;
                // mu=12.0;
                Tauion = 0.66 * 2.0 * pi * aR[im] / nlimiters / (9.79e5 * sqrt(Z / mu * Tr.Te[im]) * pow(1.0 + nevac / nr.ne[im], -0.66) * pow(1.0 + nevac / nr.nCV[im], -0.66)); //;        if (nr.ne[im]<1.0e4) {Tauion=Tauion*1.0e4/nr.ne[im];}
                dnr.dnCV[im] -= nr.nCV[im] / Tauion;
                dnr.dne[im] -= Z * nr.nCV[im] / Tauion;
                dEr.dEe[im] -= Z * nr.nCV[im] / (Tauion / gEe) * Tr.Te[im] * 3.0 / 2.0;
                dnr.dnCI[im] += nr.nCV[im] / Tauion;
            }

            if (add) ++im;
            else --im;
        }

        im = NMESHP - 1;
        add = 0;
    }
    #ifdef debug
	cout << "Limiter losses: " << "Tau = " << Tauion << " [s]" << endl;
    #endif
}

void vdrift_function() {
    double Tauion;
    double Z = 1.0;
    #ifdef debug
    	double Tauion_min = 1000.0;
        double Tauion_max = -1.0;
    #endif
//#pragma omp parallel for private(Tauion, Z)
    for (int im = 0; im < NMESHP; ++im) {
        // charge state
        // double mu=1.0; // mass unit

        Z = 1.0;                                                                                   // mu=1.0;
        Tauion = a / (100.0 * 2.0 * kb * Tr.THi[im] * 11600.0 / Z / qe / Br[im] / (aR[im] / 100)); // * sqrt(pow(Br[im],2.0)+pow(Bv,2.0)+pow(Br[im],2.0)) / (9.79e5*sqrt(Z/mu*(Tr.Te[im])) * pow(1.0+nevac/nr.ne[im],-0.66) * pow(1.0+nevac/nr.nHi[im],-0.66)) ; //;        if (nr.ne[im]<1.0e4) {Tauion=Tauion*1.0e4/nr.ne[im];}
        #ifdef debug
        	if(Tauion>Tauion_max){
			Tauion_max = Tauion;
		}
        	if(Tauion<Tauion_min){
			Tauion_min = Tauion;
		}
	#endif
        dnr.dnHi[im] -= nr.nHi[im] / Tauion;
        dEr.dEHi[im] -= Er.EHi[im] / (Tauion / gEe);
        dnr.dne[im] -= Z * nr.nHi[im] / Tauion;
        dEr.dEe[im] -= Z * nr.nHi[im] / (Tauion / gEe) * Tr.Te[im] * 3.0 / 2.0;
        // dnr.dnH2[im] += 0.5*nr.nHi[im] / Tauion;
        // dEr.dEH2[im] += 0.5*nr.nHi[im] / Tauion  *3.0/2.0*Ta0;
        Z = 1.0;                                                                                    // mu=2.0;
        Tauion = a / (100.0 * 2.0 * kb * Tr.TH2i[im] * 11600.0 / Z / qe / Br[im] / (aR[im] / 100)); // 2.0*pi*a * sqrt(pow(Br[im],2.0)+pow(Bv,2.0)+pow(Br[im],2.0)) / (9.79e5*sqrt(Z/mu*(Tr.Te[im])) * pow(1.0+nevac/nr.ne[im],-0.66) * pow(1.0+nevac/nr.nH2i[im],-0.66)) ; //;        if (nr.ne[im]<1.0e4) {Tauion=Tauion*1.0e4/nr.ne[im];}
        #ifdef debug
        	if(Tauion>Tauion_max){
			Tauion_max = Tauion;
		}
        	if(Tauion<Tauion_min){
			Tauion_min = Tauion;
		}
	#endif
        dnr.dnH2i[im] -= nr.nH2i[im] / Tauion;
        dEr.dEH2i[im] -= Er.EH2i[im] / (Tauion / gEe);
        dnr.dne[im] -= Z * nr.nH2i[im] / Tauion;
        dEr.dEe[im] -= Z * nr.nH2i[im] / (Tauion / gEe) * Tr.Te[im] * 3.0 / 2.0;
        // dnr.dnH2[im] += 1.0*nr.nH2i[im] / Tauion;
        // dEr.dEH2[im] += 1.0*nr.nH2i[im] / Tauion  *3.0/2.0*Ta0;
        Z = 1.0;                                                                                    // mu=3.0;
        Tauion = a / (100.0 * 2.0 * kb * Tr.TH3i[im] * 11600.0 / Z / qe / Br[im] / (aR[im] / 100)); // 2.0*pi*a * sqrt(pow(Br[im],2.0)+pow(Bv,2.0)+pow(Br[im],2.0)) / (9.79e5*sqrt(Z/mu*(Tr.Te[im])) * pow(1.0+nevac/nr.ne[im],-0.66) * pow(1.0+nevac/nr.nH3i[im],-0.66)) ; //;        if (nr.ne[im]<1.0e4) {Tauion=Tauion*1.0e4/nr.ne[im];}
        dnr.dnH3i[im] -= nr.nH3i[im] / Tauion;
        dEr.dEH3i[im] -= Er.EH3i[im] / (Tauion / gEe);
        dnr.dne[im] -= Z * nr.nH3i[im] / Tauion;
        dEr.dEe[im] -= Z * nr.nH3i[im] / (Tauion / gEe) * Tr.Te[im] * 3.0 / 2.0;
        //dnr.dnH2[im] += 1.5 * nr.nH3i[im] / Tauion;
        //dEr.dEH2[im] += 1.5 * nr.nH3i[im] / Tauion * 3.0 / 2.0 * Ta0;
        Z = 1.0;                                                                                     // mu=4.0;
        Tauion = a / (100.0 * 2.0 * kb * Tr.THeII[im] * 11600.0 / Z / qe / Br[im] / (aR[im] / 100)); // 2.0*pi*a * sqrt(pow(Br[im],2.0)+pow(Bv,2.0)+pow(Br[im],2.0)) / (9.79e5*sqrt(Z/mu*(Tr.Te[im])) * pow(1.0+nevac/nr.ne[im],-0.66) * pow(1.0+nevac/nr.nHeII[im],-0.66)) ; //;        if (nr.ne[im]<1.0e4) {Tauion=Tauion*1.0e4/nr.ne[im];}
        #ifdef debug
        	if(Tauion>Tauion_max){
			Tauion_max = Tauion;
		}
        	if(Tauion<Tauion_min){
			Tauion_min = Tauion;
		}
	#endif
        dnr.dnHeII[im] -= nr.nHeII[im] / Tauion;
        dEr.dEHeII[im] -= Er.EHeII[im] / (Tauion / gEe);
        dnr.dne[im] -= Z * nr.nHeII[im] / Tauion;
        dEr.dEe[im] -= Z * nr.nHeII[im] / (Tauion / gEe) * Tr.Te[im] * 3.0 / 2.0;
        // dnr.dnHeI[im] += nr.nHeII[im] / Tauion;
        // dEr.dEHeI[im] += nr.nHeII[im] / Tauion  *3.0/2.0*Ta0;
        Z = 2.0;                                                                                      // mu=4.0;
        Tauion = a / (100.0 * 2.0 * kb * Tr.THeIII[im] * 11600.0 / Z / qe / Br[im] / (aR[im] / 100)); // 2.0*pi*a * sqrt(pow(Br[im],2.0)+pow(Bv,2.0)+pow(Br[im],2.0)) / (9.79e5*sqrt(Z/mu*(Tr.Te[im])) * pow(1.0+nevac/nr.ne[im],-0.66) * pow(1.0+nevac/nr.nHeIII[im],-0.66)) ; //;        if (nr.ne[im]<1.0e4) {Tauion=Tauion*1.0e4/nr.ne[im];}
        #ifdef debug
        	if(Tauion>Tauion_max){
			Tauion_max = Tauion;
		}
        	if(Tauion<Tauion_min){
			Tauion_min = Tauion;
		}
	#endif
        dnr.dnHeIII[im] -= nr.nHeIII[im] / Tauion;
        dEr.dEHeIII[im] -= Er.EHeIII[im] / (Tauion / gEe);
        dnr.dne[im] -= Z * nr.nHeIII[im] / Tauion;
        dEr.dEe[im] -= Z * nr.nHeIII[im] / (Tauion / gEe) * Tr.Te[im] * 3.0 / 2.0;
        // dnr.dnHeI[im] += nr.nHeIII[im] / Tauion;
        // dEr.dEHeI[im] += nr.nHeIII[im] / Tauion  *3.0/2.0*Ta0;
        if (bimpur) {
            Z = 1.0;                                                                                  // mu=12.0;
            Tauion = a / (100.0 * 2.0 * kb * Tr.Te[im] * 11600.0 / Z / qe / Br[im] / (aR[im] / 100)); // 2.0*pi*a * sqrt(pow(Br[im],2.0)+pow(Bv,2.0)+pow(Br[im],2.0)) / (9.79e5*sqrt(Z/mu*(Tr.Te[im])) * pow(1.0+nevac/nr.ne[im],-0.66) * pow(1.0+nevac/nr.nCII[im],-0.66)) ; //;        if (nr.ne[im]<1.0e4) {Tauion=Tauion*1.0e4/nr.ne[im];}
            dnr.dnCII[im] -= nr.nCII[im] / Tauion;
            dnr.dne[im] -= Z * nr.nCII[im] / Tauion;
            dEr.dEe[im] -= Z * nr.nCII[im] / (Tauion / gEe) * Tr.Te[im] * 3.0 / 2.0;
            dnr.dnCI[im] += nr.nCII[im] / Tauion;
            Z = 2.0;                                                                                  // mu=12.0;
            Tauion = a / (100.0 * 2.0 * kb * Tr.Te[im] * 11600.0 / Z / qe / Br[im] / (aR[im] / 100)); // 2.0*pi*a * sqrt(pow(Br[im],2.0)+pow(Bv,2.0)+pow(Br[im],2.0)) / (9.79e5*sqrt(Z/mu*(Tr.Te[im])) * pow(1.0+nevac/nr.ne[im],-0.66) * pow(1.0+nevac/nr.nCIII[im],-0.66)) ; //;        if (nr.ne[im]<1.0e4) {Tauion=Tauion*1.0e4/nr.ne[im];}
            dnr.dnCIII[im] -= nr.nCIII[im] / Tauion;
            dnr.dne[im] -= Z * nr.nCIII[im] / Tauion;
            dEr.dEe[im] -= Z * nr.nCIII[im] / (Tauion / gEe) * Tr.Te[im] * 3.0 / 2.0;
            dnr.dnCI[im] += nr.nCIII[im] / Tauion;
            Z = 3.0;                                                                                  // mu=12.0;
            Tauion = a / (100.0 * 2.0 * kb * Tr.Te[im] * 11600.0 / Z / qe / Br[im] / (aR[im] / 100)); // 2.0*pi*a * sqrt(pow(Br[im],2.0)+pow(Bv,2.0)+pow(Br[im],2.0)) / (9.79e5*sqrt(Z/mu*(Tr.Te[im])) * pow(1.0+nevac/nr.ne[im],-0.66) * pow(1.0+nevac/nr.nCIV[im],-0.66)) ; //;        if (nr.ne[im]<1.0e4) {Tauion=Tauion*1.0e4/nr.ne[im];}
            dnr.dnCIV[im] -= nr.nCIV[im] / Tauion;
            dnr.dne[im] -= Z * nr.nCIV[im] / Tauion;
            dEr.dEe[im] -= Z * nr.nCIV[im] / (Tauion / gEe) * Tr.Te[im] * 3.0 / 2.0;
            dnr.dnCI[im] += nr.nCIV[im] / Tauion;
            Z = 4.0;                                                                                  // mu=12.0;
            Tauion = a / (100.0 * 2.0 * kb * Tr.Te[im] * 11600.0 / Z / qe / Br[im] / (aR[im] / 100)); // 2.0*pi*a * sqrt(pow(Br[im],2.0)+pow(Bv,2.0)+pow(Br[im],2.0)) / (9.79e5*sqrt(Z/mu*(Tr.Te[im])) * pow(1.0+nevac/nr.ne[im],-0.66) * pow(1.0+nevac/nr.nCV[im],-0.66)) ; //;        if (nr.ne[im]<1.0e4) {Tauion=Tauion*1.0e4/nr.ne[im];}
            dnr.dnCV[im] -= nr.nCV[im] / Tauion;
            dnr.dne[im] -= Z * nr.nCV[im] / Tauion;
            dEr.dEe[im] -= Z * nr.nCV[im] / (Tauion / gEe) * Tr.Te[im] * 3.0 / 2.0;
            dnr.dnCI[im] += nr.nCV[im] / Tauion;
        }
    }
    #ifdef debug
	cout << "Drift losses: " << "Taumin = " << Tauion_min << " [s]" << " Taumax = " << Tauion_max << endl;
    #endif
}

void bpol_function() {
    double Tauion;
    #ifdef debug
    	double Tauion_min = 1000.0;
        double Tauion_max = -1.0;
    #endif
    double Z = 1.0;  // charge state
    double Dv = 1.0; // vertical diffusion
    // double gr = 1.0; // vertical diffusion
    // double mfp = 1.0;
    double nui = 0.0;
    double mfpi = 0.0;
    double gri = 0.0;
// vertical diffusion enhanced by vertical magnetic field
//#pragma omp parallel for private(Tauion, Z, Dv)
    for (int im = 0; im < NMESHP; ++im) {
        mfpi = nr.nHi[im] * 9.79e5 * sqrt((Tr.Te[im] + Tr.THi[im] * nr.nHi[im] / nr.ne[im]) / 1.0) / max(colrate.nuHi[im], 5e3);
        gri = nr.nHi[im] * 1.02e2 / 1.0 * sqrt(1.0 * Tr.THi[im]) / Br[im] / 1e4; // cm
        nui = nr.nHi[im] * max(colrate.nuHi[im], 5e3);
        // if (im==90 & mytid==2) {cout << " 1.  " << mfpi/nr.nHi[im] << "   " << gri/nr.nHi[im] << "   " << nui << endl;}
        mfpi += nr.nH2i[im] * 9.79e5 * sqrt((Tr.Te[im] + Tr.TH2i[im] * nr.nH2i[im] / nr.ne[im]) / 2.0) / max(colrate.nuH2i[im], 5e3);
        gri += nr.nH2i[im] * 1.02e2 / 1.0 * sqrt(2.0 * Tr.TH2i[im]) / Br[im] / 1e4; // cm
        nui += nr.nH2i[im] * max(colrate.nuH2i[im], 5e3);
        // if (im==90 & mytid==2) {cout << "     " << mfpi << "   " << gri << "   " << nu << endl;}
        mfpi += nr.nH3i[im] * 9.79e5 * sqrt((Tr.Te[im] + Tr.TH3i[im] * nr.nH3i[im] / nr.ne[im]) / 3.0) / max(colrate.nuH3i[im], 5e3);
        gri += nr.nH3i[im] * 1.02e2 / 1.0 * sqrt(3.0 * Tr.TH3i[im]) / Br[im] / 1e4; // cm
        nui += nr.nH3i[im] * max(colrate.nuH3i[im], 5e3);
        // if (im==90 & mytid==2) {cout << "     " << mfpi << "   " << gri << "   " << nu << endl;}
        mfpi += nr.nHeII[im] * 9.79e5 * sqrt((Tr.Te[im] + Tr.THeII[im] * nr.nHeII[im] / nr.ne[im]) / 4.0) / max(colrate.nuHeII[im], 5e3);
        gri += nr.nHeII[im] * 1.02e2 / 1.0 * sqrt(4.0 * Tr.THeII[im]) / Br[im] / 1e4; // cm
        nui += nr.nHeII[im] * max(colrate.nuHeII[im], 5e3);
        // if (im==90 & mytid==2) {cout << "     " << mfpi << "   " << gri << "   " << nu << endl;}
        mfpi += nr.nHeIII[im] * 9.79e5 * sqrt(2.0 * (Tr.Te[im] + Tr.THeIII[im] * nr.nHeIII[im] / nr.ne[im]) / 4.0) / max(colrate.nuHeIII[im], 5e3);
        gri += nr.nHeIII[im] * 1.02e2 / 2.0 * sqrt(4.0 * Tr.THeIII[im]) / Br[im] / 1e4; // cm
        nui += nr.nHeIII[im] * max(colrate.nuHeIII[im], 5e3);
        // if (im==90 & mytid==2) {cout << "     " << mfpi << "   " << gri << "   " << nu << endl;}
        mfpi = mfpi / (nr.nHi[im] + nr.nH2i[im] + nr.nH3i[im] + nr.nHeII[im] + nr.nHeIII[im]);
        gri = gri / (nr.nHi[im] + nr.nH2i[im] + nr.nH3i[im] + nr.nHeII[im] + nr.nHeIII[im]);
        nui = nui / (nr.nHi[im] + nr.nH2i[im] + nr.nH3i[im] + nr.nHeII[im] + nr.nHeIII[im]);

        Dv = Dfsave * 0.333 * nui * mfpi * (gri + mfpi * Bv / Br[im]);
        // Tauion = pow(min(b,-1.1696*aR[im]+187.8655),2.0)/(2.0*Dv);
        Tauion = pow(b, 2.0) / (2.0 * Dv);
        #ifdef debug
        	if(Tauion>Tauion_max){
			Tauion_max = Tauion;
		}
        	if(Tauion<Tauion_min){
			Tauion_min = Tauion;
		}
	#endif

        Z = 1.0; // mu=1.0;
        // gr = 1.02e2*sqrt(mu*Tr.THi[im])/Br[im]/1e4;
        // mfp= 1.0/.colrate.nuHi[im] * 9.79e5*sqrt(Z/mu*(Tr.Te[im]+Tr.THi[im]*nr.nHi[im]/nr.ne[im])) * pow(1.0+nevac/nr.ne[im],-0.66) * pow(1.0+nevac/nr.nHi[im],-0.66) ;
        // Dv = Dfsave * 0.333 * .colrate.nuHi[im] * mfp * (gr + mfp*Bv/Br[im] );
        // Tauion = pow(min(b,-1.1696*aR[im]+187.8655),2.0)/(2.0*Dv);
        dnr.dnHi[im] -= nr.nHi[im] / Tauion;
        dEr.dEHi[im] -= Er.EHi[im] / (Tauion / gEd);
        dnr.dne[im] -= Z * nr.nHi[im] / Tauion;
        dEr.dEe[im] -= Z * nr.nHi[im] / (Tauion / gEd) * Tr.Te[im] * 3.0 / 2.0;
        // dnr.dnH2[im] += 0.5*nr.nHi[im] / Tauion;
        // dEr.dEH2[im] += 0.5*nr.nHi[im] / Tauion  *3.0/2.0*Ta0;
        Z = 1.0; // mu=2.0;
        // gr = 1.02e2*sqrt(mu*Tr.TH2i[im])/Br[im]/1e4;
        // mfp= 1.0/.colrate.nuH2i[im] * 9.79e5*sqrt(Z/mu*(Tr.Te[im]+Tr.TH2i[im]*nr.nH2i[im]/nr.ne[im])) * pow(1.0+nevac/nr.ne[im],-0.66) * pow(1.0+nevac/nr.nH2i[im],-0.66) ;
        // Dv = Dfsave * 0.333 * .colrate.nuH2i[im] * mfp * (gr + mfp*Bv/Br[im] );
        // Tauion = pow(min(b,-1.1696*aR[im]+187.8655),2.0)/(2.0*Dv);
        dnr.dnH2i[im] -= nr.nH2i[im] / Tauion;
        dEr.dEH2i[im] -= Er.EH2i[im] / (Tauion / gEd);
        dnr.dne[im] -= Z * nr.nH2i[im] / Tauion;
        dEr.dEe[im] -= Z * nr.nH2i[im] / (Tauion / gEd) * Tr.Te[im] * 3.0 / 2.0;
        // dnr.dnH2[im] += 1.0*nr.nH2i[im] / Tauion;
        // dEr.dEH2[im] += 1.0*nr.nH2i[im] / Tauion  *3.0/2.0*Ta0;
        Z = 1.0;
        // mu=3.0;
        //  gr = 1.02e2*sqrt(mu*Tr.TH3i[im])/Br[im]/1e4;
        //  mfp= 1.0/.colrate.nuH3i[im] * 9.79e5*sqrt(Z/mu*(Tr.Te[im]+Tr.TH3i[im]*nr.nH3i[im]/nr.ne[im])) * pow(1.0+nevac/nr.ne[im],-0.66) * pow(1.0+nevac/nr.nH3i[im],-0.66) ;
        //  Dv = Dfsave * 0.333 * .colrate.nuH3i[im] * mfp * (gr + mfp*Bv/Br[im] );
        //  Tauion = pow(min(b,-1.1696*aR[im]+187.8655),2.0)/(2.0*Dv);
        dnr.dnH3i[im] -= nr.nH3i[im] / Tauion;
        dEr.dEH3i[im] -= Er.EH3i[im] / (Tauion / gEd);
        dnr.dne[im] -= Z * nr.nH3i[im] / Tauion;
        dEr.dEe[im] -= Z * nr.nH3i[im] / (Tauion / gEd) * Tr.Te[im] * 3.0 / 2.0;
        // dnr.dnH2[im] += 1.5*nr.nH3i[im] / Tauion;
        // dEr.dEH2[im] += 1.5*nr.nH3i[im] / Tauion  *3.0/2.0*Ta0;
        Z = 1.0;
        // mu=4.0;
        //  gr = 1.02e2*sqrt(mu*Tr.THeII[im])/Br[im]/1e4;
        //  mfp= 1.0/.colrate.nuHeII[im] * 9.79e5*sqrt(Z/mu*(Tr.Te[im]+Tr.THeII[im]*nr.nHeII[im]/nr.ne[im])) * pow(1.0+nevac/nr.ne[im],-0.66) * pow(1.0+nevac/nr.nHeII[im],-0.66) ;
        //  Dv = Dfsave * 0.333 * .colrate.nuHeII[im] * mfp * (gr + mfp*Bv/Br[im] );
        //  Tauion = pow(min(b,-1.1696*aR[im]+187.8655),2.0)/(2.0*Dv);
        //  if (im==90 ) {cout << "  " <<  mytid << "  " << .colrate.nuHi[im] << "  " << Tauion << "  " << Dv  << "  " << Tr.THi[im] << endl;}
        dnr.dnHeII[im] -= nr.nHeII[im] / Tauion;
        dEr.dEHeII[im] -= Er.EHeII[im] / (Tauion / gEd);
        dnr.dne[im] -= Z * nr.nHeII[im] / Tauion;
        dEr.dEe[im] -= Z * nr.nHeII[im] / (Tauion / gEd) * Tr.Te[im] * 3.0 / 2.0;
        // dnr.dnHeI[im] += nr.nHeII[im] / Tauion;
        // dEr.dEHeI[im] += nr.nHeII[im] / Tauion  *3.0/2.0*Ta0;
        Z = 2.0;
        // gr = 1.02e2*sqrt(mu*Tr.THeIII[im])/Br[im]/1e4;
        // mfp= 1.0/.colrate.nuHeIII[im] * 9.79e5*sqrt(Z/mu*(Tr.Te[im]+Tr.THeIII[im]*nr.nHeIII[im]/nr.ne[im])) * pow(1.0+nevac/nr.ne[im],-0.66) * pow(1.0+nevac/nr.nHeIII[im],-0.66) ;
        // Dv = Dfsave * 0.333 * .colrate.nuHeIII[im] * mfp * (gr + mfp*Bv/Br[im] );
        // Tauion = pow(min(b,-1.1696*aR[im]+187.8655),2.0)/(2.0*Dv);
        dnr.dnHeIII[im] -= nr.nHeIII[im] / Tauion;
        dEr.dEHeIII[im] -= Er.EHeIII[im] / (Tauion / gEd);
        dnr.dne[im] -= Z * nr.nHeIII[im] / Tauion;
        dEr.dEe[im] -= Z * nr.nHeIII[im] / (Tauion / gEd) * Tr.Te[im] * 3.0 / 2.0;
        // dnr.dnHeI[im] += nr.nHeIII[im] / Tauion;
        // dEr.dEHeI[im] += nr.nHeIII[im] / Tauion  *3.0/2.0*Ta0;
        if (bimpur) {
            Z = 1.0; // mu=12.0;
            // gr = 1.02e2*sqrt(mu*Tr.Te[im])/Br[im]/1e4;
            // mfp= 1.0/5e4 * 9.79e5*sqrt(Z/mu*(Tr.Te[im])) * pow(1.0+nevac/nr.ne[im],-0.66) * pow(1.0+nevac/nr.nCII[im],-0.66) ;
            // Dv = Dfsave * 0.333 * 5e4 * mfp * (gr + mfp*Bv/Br[im] );
            // Tauion = pow(min(b,-1.1696*aR[im]+187.8655),2.0)/(2.0*Dv);
            dnr.dnCII[im] -= nr.nCII[im] / Tauion;
            dnr.dne[im] -= Z * nr.nCII[im] / Tauion;
            dEr.dEe[im] -= Z * nr.nCII[im] / (Tauion / gEd) * Tr.Te[im] * 3.0 / 2.0;
            dnr.dnCI[im] += nr.nCII[im] / Tauion;
            Z = 2.0;
            // mu=12.0;
            //  gr = 1.02e2*sqrt(mu*Tr.Te[im])/Br[im]/1e4;
            //  mfp= 1.0/5e4 * 9.79e5*sqrt(Z/mu*(Tr.Te[im])) * pow(1.0+nevac/nr.ne[im],-0.66) * pow(1.0+nevac/nr.nCIII[im],-0.66) ;
            //  Dv = Dfsave * 0.333 * 5e4 * mfp * (gr + mfp*Bv/Br[im] );
            //  Tauion = pow(min(b,-1.1696*aR[im]+187.8655),2.0)/(2.0*Dv);
            dnr.dnCIII[im] -= nr.nCIII[im] / Tauion;
            dnr.dne[im] -= Z * nr.nCIII[im] / Tauion;
            dEr.dEe[im] -= Z * nr.nCIII[im] / (Tauion / gEd) * Tr.Te[im] * 3.0 / 2.0;
            dnr.dnCI[im] += nr.nCIII[im] / Tauion;
            Z = 3.0;
            // mu=12.0;
            //  gr = 1.02e2*sqrt(mu*Tr.Te[im])/Br[im]/1e4;
            //  mfp= 1.0/5e4 * 9.79e5*sqrt(Z/mu*(Tr.Te[im])) * pow(1.0+nevac/nr.ne[im],-0.66) * pow(1.0+nevac/nr.nCIV[im],-0.66) ;
            //  Dv = Dfsave * 0.333 * 5e4 * mfp * (gr + mfp*Bv/Br[im] );
            //  Tauion = pow(min(b,-1.1696*aR[im]+187.8655),2.0)/(2.0*Dv);
            dnr.dnCIV[im] -= nr.nCIV[im] / Tauion;
            dnr.dne[im] -= Z * nr.nCIV[im] / Tauion;
            dEr.dEe[im] -= Z * nr.nCIV[im] / (Tauion / gEd) * Tr.Te[im] * 3.0 / 2.0;
            dnr.dnCI[im] += nr.nCIV[im] / Tauion;
            Z = 4.0;
            // mu=12.0;
            //  gr = 1.02e2*sqrt(mu*Tr.Te[im])/Br[im]/1e4;
            //  mfp= 1.0/5e4 * 9.79e5*sqrt(Z/mu*(Tr.Te[im])) * pow(1.0+nevac/nr.ne[im],-0.66) * pow(1.0+nevac/nr.nCV[im],-0.66) ;
            //  Dv = Dfsave * 0.333 * 5e4 * mfp * (gr + mfp*Bv/Br[im] );
            //  Tauion = pow(min(b,-1.1696*aR[im]+187.8655),2.0)/(2.0*Dv);
            dnr.dnCV[im] -= nr.nCV[im] / Tauion;
            dnr.dne[im] -= Z * nr.nCV[im] / Tauion;
            dEr.dEe[im] -= Z * nr.nCV[im] / (Tauion / gEd) * Tr.Te[im] * 3.0 / 2.0;
            dnr.dnCI[im] += nr.nCV[im] / Tauion;
        }
    }

    #ifdef debug
	cout << "Poloidal losses: " << "Dv = " << Dv/1e4 << " [m2/s]; Taumin = " << Tauion_min << " [s]; Taumax = " << Tauion_max << endl;
    #endif
}
