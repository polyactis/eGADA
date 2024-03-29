/*=================================================================
 * BaseGADA.h
 * Header file BaseGADA.c.
 *=================================================================*/
/*
 This File is part of eGADA

 eGADA: enhanced Genome Alteration Detection Algorithm

 eGADA is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 any later version.

 eGADA is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with eGADA.  If not, see <http://www.gnu.org/licenses/>.

 Authors:
 Yu S. Huang  polyactis@gmail.com
 Roger Pique-Regi    piquereg@usc.edu

 */

#ifndef _BaseGADA_H_
#define _BaseGADA_H_

//#include "matlabdefines.h"
#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <map>	//for hash_map
#include <set>	//for set
// to customize hash
#include <functional>
// yh: to customize boost::hash
#include <boost/functional/hash.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
// red-black tree to store segment breakpoint, score, etc.
#include "RedBlackTree.h"


#ifndef GADABIN
#include <boost/python.hpp>
#endif
using namespace std;


class BreakPointKey{
    /*
     * class defines the key of a breakpoint in a red-black tree.
     */

    // for output formatting
    friend ostream& operator<<(ostream& out, const BreakPointKey& bpKey){
        /*
         * something wrong here, it can't be streamed to an ostream
         */
        out << boost::format("position=%1%, tscore=%2%, weight=%3%, length=%4%, "
            "MinSegLen=%5%, totalLength=%6%")%
                bpKey.position % bpKey.tscore % bpKey.weight %
                bpKey.segmentLength % bpKey.MinSegLen % bpKey.totalLength;
        return out;
    }
    /* The comparison operator for RBTree nodes determines how to order them.
     * If at least one t-stat/t-score of two breakpoints is below the threshold,
     *   order them based on t-score.
     * Otherwise, order them by segment length, which is defined as the length
     *   of the shorter flanking fragment.
    */
    friend bool operator<(const BreakPointKey& a, const BreakPointKey& b){
        /*
         * make sure no duplicate (or close keys, for floating values) keys 
         *      exist in red-black trees.
         */

        //both below minTScore, then order them by t-score
        if(a.tscore<a.minTScore && b.tscore<b.minTScore){
            if (a.tscore==b.tscore){
                return a.segmentLength<b.segmentLength;
            }
            else{
                return a.tscore<b.tscore;
            }
        }
        //one below, one above
        else if ((a.tscore<a.minTScore && b.tscore>=b.minTScore) || 
            (a.tscore>=a.minTScore && b.tscore<b.minTScore)){
            if (a.tscore==b.tscore){
                return a.segmentLength<b.segmentLength;
            }
            else{
                //order by tscore
                return a.tscore<b.tscore;
            }
        }
        //both above minTScore or one equal one above or both equal), 
        else{
            if (a.segmentLength==b.segmentLength){
                return a.tscore<b.tscore;
            }
            else{
                return a.segmentLength<b.segmentLength;
            }
        }

    }

    friend bool operator>(const BreakPointKey& a, const BreakPointKey& b){
        /*
         * make sure no duplicate (or close keys, for floating values) keys 
         *  exist in red-black trees.
         */

        //both below minTScore, then order them by score
        if(a.tscore<a.minTScore && b.tscore<b.minTScore){
            if (a.tscore==b.tscore){
                return a.segmentLength>b.segmentLength;
            }
            else{
                return a.tscore>b.tscore;
            }
        }
        //one below, one above
        else if ((a.tscore<a.minTScore && b.tscore>=b.minTScore) || 
            (a.tscore>=a.minTScore && b.tscore<b.minTScore)){
            if (a.tscore==b.tscore){
                return a.segmentLength>b.segmentLength;
            }
            else{
                //order by tscore
                return a.tscore>b.tscore;
            }
        }
        //both above minTScore or one equal one above or both equal), 
        else{
            if (a.segmentLength==b.segmentLength){
                return a.tscore>b.tscore;
            }
            else{
                return a.segmentLength>b.segmentLength;
            }
        }
    }

    friend bool operator==(const BreakPointKey& a, const BreakPointKey& b) {
        return a.tscore == b.tscore;
    }

public:
    long position;
    double weight;
    double tscore;
    // length of the shorter flanking segment.
    long segmentLength;
    long MinSegLen;
    double minTScore;
    long totalLength;
    void* nodePtr;
    
    BreakPointKey(){
        position=1;
        weight=0;
        tscore=0;
        segmentLength=0;
        MinSegLen=0;
        minTScore = 0;
        totalLength=0;
        nodePtr = NULL;
    }
    
    BreakPointKey(long _position, double _weight, double _tscore, long _segmentLength, 
        long _MinSegLen, double _minTScore, long _totalLength):
        position(_position),weight(_weight), tscore(_tscore), 
        segmentLength(_segmentLength), MinSegLen(_MinSegLen), 
        minTScore(_minTScore), totalLength(_totalLength){
        nodePtr = NULL;
    };

    ~BreakPointKey(){
    }
};


class BreakPoint{
    // for output
    friend ostream& operator<<(ostream& out, const BreakPoint& breakPoint){
        out << boost::format("position=%1%, tscore=%2%, weight=%3%, length=%4%, "
            "MinSegLen=%5%, totalLength=%6%")%
                breakPoint.position % breakPoint.tscore % breakPoint.weight %
                breakPoint.segmentLength % breakPoint.MinSegLen % breakPoint.totalLength;
        return out;
    }
    friend bool operator<(const BreakPoint& a, const BreakPoint& b){
        /*
         * to be used in std::sort
         */
        return a.position<b.position;
    }

    friend bool operator==(const BreakPoint& a, const BreakPoint& b) {
        return a.position==b.position;
    }

public:
    long position;
    double weight;
    double tscore;
    long segmentLength;
    long MinSegLen;
    double minTScore;
    long totalLength;
    BreakPoint* leftBreakPointPtr;
    BreakPoint* rightBreakPointPtr;
    void* nodePtr;
    BreakPoint(){
        position=1;
        weight=0;
        tscore=0;
        segmentLength=0;
        MinSegLen=0;
        totalLength=0;
        leftBreakPointPtr=NULL;
        rightBreakPointPtr=NULL;
        nodePtr = NULL;
    }
    BreakPoint(long _position, double _weight, double _tscore, long _segmentLength, 
            long _MinSegLen, double _minTScore, long _totalLength):
        position(_position),weight(_weight), tscore(_tscore), 
        segmentLength(_segmentLength), MinSegLen(_MinSegLen), minTScore(_minTScore), 
        totalLength(_totalLength){
        leftBreakPointPtr=NULL;
        rightBreakPointPtr=NULL;
        nodePtr = NULL;
    };

    ~BreakPoint(){};
    
    void setLeftBreakPoint(BreakPoint* bpPtr){
        this->leftBreakPointPtr=bpPtr;
    }
    
    BreakPoint* getLeftBreakPoint(){
        return this->leftBreakPointPtr;
    }
    
    void setRightBreakPoint(BreakPoint* bpPtr){
        this->rightBreakPointPtr=bpPtr;
    }
    
    BreakPoint* getRightBreakPoint(){
        return this->rightBreakPointPtr;
    }
    
    void setRBTreeNodePtr(void* ndPtr){
        this->nodePtr = ndPtr;
    }
    
    void* getRBTreeNodePtr(){
        return this->nodePtr;
    }
    
    BreakPointKey* getKeyPointer(){
        BreakPointKey* bpKeyPtr = new BreakPointKey(position, weight, tscore, 
            segmentLength, MinSegLen, minTScore, totalLength);
        return bpKeyPtr;
    }
    
    BreakPointKey getKey(){
        BreakPointKey bpKey = BreakPointKey(position, weight, tscore, 
            segmentLength, MinSegLen, minTScore, totalLength);
        return bpKey;
    }
    
    void removeItself(){
        /*
         * updating left and right breakpoint
         */
        long j;
        double iC, iL, iR, M;
        double h0;

        M = (double) totalLength;
        //initialize
        BreakPoint* leftLeftBreakPointPtr = NULL;
        BreakPoint* rightRightBreakPointPtr = NULL;

        if (leftBreakPointPtr!=NULL && rightBreakPointPtr!=NULL){
            iL = (double) leftBreakPointPtr->position;
            iC = (double) this->position;
            iR = (double) rightBreakPointPtr->position;

            leftLeftBreakPointPtr = leftBreakPointPtr->leftBreakPointPtr;
            rightRightBreakPointPtr = rightBreakPointPtr->rightBreakPointPtr;
            //update the weights first
            
            //leftBreakPointPtr is NOT the left most.
            if (leftLeftBreakPointPtr!=NULL){
                leftBreakPointPtr->weight = leftBreakPointPtr->weight
                    + sqrt((M - iL) / (M - iC) * iL / iC) * (iR - iC) / (iR - iL) *weight;
            }
            //rightBreakPointPtr is NOT the right most break point
            if (rightRightBreakPointPtr!=NULL){
                rightBreakPointPtr->weight = rightBreakPointPtr->weight
                    + sqrt((M - iR) / (M - iC) * iR / iC) * (iC - iL) / (iR - iL) * weight;
            }
            //update the tscoe and segmentLength, which needs the updated weights
            if (leftLeftBreakPointPtr!=NULL){
                h0 = (double) (M - leftBreakPointPtr->position) * 
                    (double) leftBreakPointPtr->position / M * 
                    (double) (rightBreakPointPtr->position - leftLeftBreakPointPtr->position) / 
                    (double) (rightBreakPointPtr->position - leftBreakPointPtr->position) / 
                    (double) (leftBreakPointPtr->position - leftLeftBreakPointPtr->position);
                
                leftBreakPointPtr->tscore = fabs(leftBreakPointPtr->weight) / sqrt(h0);
                leftBreakPointPtr->segmentLength = min(leftBreakPointPtr->position - 
                        leftLeftBreakPointPtr->position ,
                    rightBreakPointPtr->position-leftBreakPointPtr->position);
            }
            if (rightRightBreakPointPtr!=NULL){
                h0 = (double) (M - rightBreakPointPtr->position) * 
                    (double) rightBreakPointPtr->position / M * 
                    (double) (rightRightBreakPointPtr->position - leftBreakPointPtr->position) / 
                    (double) (rightRightBreakPointPtr->position - rightBreakPointPtr->position) / 
                    (double) (rightBreakPointPtr->position - leftBreakPointPtr->position);
                
                rightBreakPointPtr->tscore = fabs(rightBreakPointPtr->weight) / sqrt(h0);
                rightBreakPointPtr->segmentLength= min(
                    rightRightBreakPointPtr->position - rightBreakPointPtr->position,
                    rightBreakPointPtr->position-leftBreakPointPtr->position);
            }
        }
        //update the left & right of the left & right
        if (leftBreakPointPtr!=NULL){
            leftBreakPointPtr->rightBreakPointPtr = rightBreakPointPtr ;
        }
        if (rightBreakPointPtr!=NULL){
            rightBreakPointPtr->leftBreakPointPtr = leftBreakPointPtr;
        }
        //release the memory
        //delete leftBreakPointPtr;
        //delete rightBreakPointPtr;
        //delete nodePtr;
    }
    //define methods to compute /update score/length/weight/position 
    //  based on neighboring breakpoints
};


namespace std {
    template<> struct hash<BreakPoint> {
        /*
         * the hash function for BreakPoint, this requires g++ flag "-std=c++0x"
         */
        size_t operator()(const BreakPoint& bp) const {

            // Start with a hash value of 0    .
            size_t seed = 0;

            // Modify 'seed' by XORing and bit-shifting in
            // one member of 'Key' after the other:
            boost::hash_combine(seed, boost::hash_value(bp.position));
            //boost::hash_combine(seed, boost::hash_value(bp.segmentLength));
            return seed;
        }
    };

};

typedef set<BreakPoint*> rbNodeDataType;
typedef RedBlackTree<BreakPointKey, rbNodeDataType > treeType;
typedef RedBlackTreeNode<BreakPointKey, rbNodeDataType > rbNodeType;

class BaseGADA{


public:
    double *Wext; //IO Breakpoint weights extended notation...
    long *Iext; //IO Breakpoint positions in extended notation...
    double *tscore;
    long *pK; //IO Number breakpoint positions remaining.
    double T; //IP  Threshold to prune
    long MinSegLen;	//minimum segment length

    //verbosity... set equal to 1 to see messages of SBLandBE(). 0 to not see them
    long debug;
    int report;
    // total/max length
    long M;

    long i;
    long K;
    //would store normalized array of input data (raw-mean)
    double *tn;
    double *inputDataArray;
    long *SegLen;
    double *SegAmp;
    double *SegState;
    double delta;
    long numEMsteps;
    long noOfBreakpointsAfterSBL;

    double *alpha, *aux;

    //double T2=5.0; //Segment collapse to base Non-Alteration level threshold
    double BaseAmp; //Base-level
    double a; //SBL hyperprior parameter
    //Variance observed, if negative value, it will be estimated by the mean of the differences
    // I would recommend to be estimated on all the chromosomes and as a trimmed mean.
    double sigma2;
    
    long SelectClassifySegments; //Classify segment into altered state (1), otherwise 0
    long SelectEstimateBaseAmp; //Estimate Neutral hybridization amplitude.

    long maxNoOfIterations;	//=50000, //10000 is enough usually
    //1E-10 or 1E-8 seems to work well for this parameter. -- => ++ conv time
    //1E8 better than 1E10 seems to work well for this parameter. -- => -- conv time
    double convergenceDelta;
    double convergenceMaxAlpha;	// 1E8 Maximum number of iterations to reach convergence...
    double convergenceB;	// a number related to convergence = 1E-20
    double ymean;	//mean of inputDataArray
    // how often to report progress during backward elimination, default is 100K
    int reportIntervalDuringBE;

    BaseGADA(double* _inputDataArray, long _M, double _sigma2, double _BaseAmp, 
        double _a, double _T, long _MinSegLen,
        long _debug , double _convergenceDelta,
        long _maxNoOfIterations, double _convergenceMaxAlpha, double _convergenceB, 
        int _reportIntervalDuringBE):
        
        inputDataArray(_inputDataArray), M(_M),  sigma2(_sigma2), BaseAmp(_BaseAmp), 
        a(_a), T(_T), MinSegLen(_MinSegLen),
        debug(_debug),
        convergenceDelta(_convergenceDelta),
        maxNoOfIterations (_maxNoOfIterations), convergenceMaxAlpha(_convergenceMaxAlpha),
        convergenceB(_convergenceB), reportIntervalDuringBE(_reportIntervalDuringBE){
        noOfBreakpointsAfterSBL = 0;
    }

    ~BaseGADA(){
#ifndef GADALib 
//this causes eGADA::run() of eGADA.so to crash right before returning segments to python
        free(SegLen);
        free(SegAmp);
        free(SegState);
#endif
    }

    void reconstruct(double *wr, long M, double *aux_vec);
    void BubbleSort(long *I, long L);
    void doubleBubbleSort(double *D, long *I, long L);
    void TrisolveREG(double *t0, double *tu, double *tl, double *coef, double *sol,
            long sizeh0);
    void DiagOfTriXTri(double *ll, double *l0, double *lu, double *rl, double *r0,
            double *ru, double *d, long N);
    void tridiagofinverse(double *t0, double *tl, double *itl, double *it0,
            double *itu, long N, double *d, double *e);
    void ForwardElimination(double *A, long N);
    void BackSubstitution(double *A, long N);
    void BackwardElimination(double *A, long N);
    void TriSolveINV(double *AA, long M, long N, double *x, double *d, double *e);
    void ComputeH(double *h0, double *h1, long M);
    void ComputeFdualXb(long M, double *convergenceB);
    // REMOVED void ComputeHs(long *s,double *a,long M,long Ms,double *h0,double *h1);
    void ComputeHs(long *s, long M, long Ms, double *h0, double *h1);
    void TriSymGaxpy(double *t0, double *t1, double *x, long M, double *y);
    void ComputeT(double *h0, double *h1, long M, double *alfa, double sigma,
            double *t0, double *tl, double *tu);
    long findminus(double *alpha, long Ms, double convergenceMaxAlpha, long *sel);
    long simpletresholding(double *inputvector, long N, double thres, double *disc);
    void computesegmentmeans(double *inputvector, long N, double *disc,
            long numdisc, double *amp);
    void reconstructoutput(double *rec, long N, double *disc, long numdisc,
            double *amp);
    
    long SBL(double *y, //I -- 1D array with the input signal
            long *I, //IO -- 1D array with the initial (final) candidate breakpoints
            double *alpha, //I -- 1D array with the initial (final) hyperparameter inv. varainces.
            double *w, //O -- 1D array containing the breakpoint weigths or posterior mean.
            double *sigw, //O -- Posterior variances, I would not really need them since they can be computed from alpha and H
            long M, //Initial size of the array in y
            long *K, //Size of the I alpha w

            //Algorithm parameters:
            double sigma2, //Noise estimated
            double a, //
            double convergenceB, double convergenceMaxAlpha, //Basis reduction parameter
            long maxNoOfIterations, //Max number of iterations
            double convergenceDelta, //Tolerance for convergence
            long debug //verbosity... set equal to 1 to see messages  0 to not see them
            );

    //To eliminate...
    long BEthresh(
            double *Scores, long Nscores, double *wr, long *indsel,
            long *pointNumRem, double *pointTau);
    
    //Returns breakpoint list lenght.
    long SBLandBE();

    void Project(double *y, long M, long *I, long L, double *xI, double *wI);
    void IextToSegLen();
    void IextWextToSegAmp();

    // computes z=F'y for entire possilbe breakpoint positions (normalized PWC)
    void CompZ(double *y, double *z, long M);
    
    void ComputeHsIext(
    //input variables:
            long *Iext, // Indices of selection,
            long K, // Length of the indices,
            double *h0, // Returning diagonal of H,
            double *h1 // Returning upper diagonal of H
            );
    
    void ProjectCoeff( //IextYobs2Wext
            double *y, long M, long *Iext, long K, double *Wext);
    
    //Uses a T test to decide which segments collapse to neutral
    void CollapseAmpTtest();

    // Returns BaseAmp corresponding to the base level.
    double CompBaseAmpMedianMethod();

    void ClassifySegments(double *SegAmp, long *SegLen, double *SegState, long K,
            double BaseAmp, double ploidy, double sigma2, //Reference noise
            double T //Critical value that decides when to colapse
            );

    void ComputeTScores(const double *Wext, const long *Iext, double *Scores, long K,
            long start, long end);

    long BEwTscore(double *Wext, //IO Breakpoint weights extended notation...
            long *Iext, //IO Breakpoint positions in extended notation...
            double *tscore, long *pK, //IO Number breakpoint positions remaining.
            double T, //IP  Threshold to prune
            long MinSegLen=0,	//minimum segment length
            long debug=0
            );

    //Returns breakpoint list lenght. with T and MinSegLen
    long BEwTandMinLen(
            double *Wext, //IO Breakpoint weights extended notation...
            long *Iext, //IO Breakpoint positions in extended notation...
            long *pK, //IO Number breakpoint positions remaining.
            double sigma2, //IP If sigma2
            double T, //IP  Threshold to prune,  T=T*sqrt(sigma2);
            long MinSegLen, //IP Minimum length of the segment.
            long debug
            );
    long RemoveBreakpoint(double *Wext, long *Iext, double *tscore, long K, 
        long indexOfSegmentToRemove);

};

#endif //_BaseGADA_H_
