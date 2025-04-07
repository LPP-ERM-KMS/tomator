#ifndef REACTIONRATES_H
#define REACTIONRATES_H

#include "constants.h"
#include <cmath>
#include <stdio.h>
#include <unordered_map>
#include <string.h>

// Reaction labels
enum Reaction {
    REAC211 = 211,
    REAC212 = 212,
    REAC213 = 213,
    REAC214a = 21401,
    REAC214b = 21402,
    REAC215 = 215,
    REAC216 = 216,
    REAC217 = 217,
    REAC218a = 21801,
    REAC218b = 21802,
    REAC221a = 22101,
    REAC221b = 22102,
    REAC222 = 222,
    REAC223 = 223,
    REAC224 = 224,
    REAC225 = 225,
    REAC226 = 226,
    REAC227 = 227,
    REAC228 = 228,
    REAC229 = 229,
    REAC2210 = 2210,
    REAC2211 = 2211,
    REAC2212 = 2212,
    REAC2213 = 2213,
    REAC2214 = 2214,
    REAC2215 = 2215,
    REAC2216 = 2216,
    REAC2217 = 2217,
    DISS = 1,
    IONI = 2,
    RECO = 3,
    ELAS = 4,
    IHE1 = 10,
    IHE2 = 11,
    RHE2 = 12,
    RHE3 = 13,
    REAC311 = 311,
    REAC312 = 312,
    REAC316 = 316,
    REAC321 = 321,
    REAC322 = 322,
    REAC323 = 323,
    REAC325 = 325,
    REAC326 = 326,
    REAC332 = 332,
    REAC431 = 431,
    REAC433 = 433,
    REAC441 = 441,
    REAC521 = 521,
    REAC523 = 523,
    REAC531 = 531,
    REAC532 = 532,
    REAC621 = 621,
    REAC631 = 631,
    CXHe3H = 20,
    CXHe2H = 21,
    CXHe3He1 = 22,
    HH2 = 30,
    HiH = 31,
    H2iH = 32,
    H3iH = 33,
    HeIIH = 34,
    HiH2 = 35,
    H2iH2 = 36,
    H3iH2 = 37,
    HeIIH2 = 38,
    HHe = 39,
    HiHe = 40,
    HeIIHe = 41,
    HeIHeI = 42,
    H2H2 = 43,
    HH = 44,
    IC1 = 50,
    IC2 = 51,
    IC3 = 52,
    IC4 = 53
};

struct ReactionType {
    int chapter;
    int section;
    int subsection;
    char type; // '\0' if no type, otherwise 'a', 'b', etc.
};

struct ReactionData {
    double values[9]; // Assuming 9 values per reaction
};

//////////////////////////////////////////////////////////////////////////////////
/// Help functions
//////////////////////////////////////////////////////////////////////////////////
void extractData(ReactionData *reaction, FILE *texFile);
void initReactionDataMap(char *fileName);
double kcalc(double Rfit[9], double T);
double kcalc2d(double Rfit2d[9][9], double E, double T);
int findindex(double arr[], int size, double val);
double interpolate1D(double arrx[], double arrz[], int sizex, double valx);
template <int sizex, int sizey>
double interpolate2D(double arrx[sizex], double arry[sizey], double arrz[sizey][sizex], double valx, double valy) {
    // 2D linear interpolation
    double zip0, zip1, zip;
    int indx, indy;
    indx = 0;
    indy = 0;
    indx = findindex(arrx, sizex, valx);
    indy = findindex(arry, sizey, valy);
    if (valx <= arrx[0]) {
        valx = arrx[0];
        indx = 0;
    } else if (valx >= arrx[sizex - 1]) {
        valx = arrx[sizex - 1];
        indx = sizex - 2;
    }
    if (valy <= arry[0]) {
        valy = arry[0];
        indy = 0;
    } else if (valy >= arry[sizey - 1]) {
        valy = arry[sizey - 1];
        indy = sizey - 2;
    }

    zip0 = arrz[indy][indx] + (arrz[indy][indx + 1] - arrz[indy][indx]) * (valx - arrx[indx]) / (arrx[indx + 1] - arrx[indx]);
    zip1 = arrz[indy + 1][indx] + (arrz[indy + 1][indx + 1] - arrz[indy + 1][indx]) * (valx - arrx[indx]) / (arrx[indx + 1] - arrx[indx]);
    zip = zip0 + (zip1 - zip0) * (valy - arry[indy]) / (arry[indy + 1] - arry[indy]);

    return zip;
}

//////////////////////////////////////////////////////////////////////////////////
/// Collisions based on REITER/JANVEV
//////////////////////////////////////////////////////////////////////////////////
double RR(int reac, double E, double T);

//////////////////////////////////////////////////////////////////////////////////
/// Collisions based on DIRK W.'s calculations
//////////////////////////////////////////////////////////////////////////////////
double RRH2(int reac, double ne, double Te);

//////////////////////////////////////////////////////////////////////////////////
/// Collisions with He based on SASCHA M's calculations
//////////////////////////////////////////////////////////////////////////////////
double RRHe(int reac, double ne, double Te);
double RRCX(int reac, double T1, double T2);

//////////////////////////////////////////////////////////////////////////////////
/// Ion collisions based on own calculations with REITER/JAVNEV cross sections
//////////////////////////////////////////////////////////////////////////////////
double RRion(int reac, double T1, double T2);

//////////////////////////////////////////////////////////////////////////////////
/// Elastic collisions
//////////////////////////////////////////////////////////////////////////////////
double RRel(int reac, double T1, double T2);

//////////////////////////////////////////////////////////////////////////////////
/// Carbon
//////////////////////////////////////////////////////////////////////////////////
double RRC(int reac, double T);

#endif // REACTIONRATES_h