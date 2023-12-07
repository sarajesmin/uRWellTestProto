/*
 * File:   uRwellTools.h
 * Author: rafopar
 *
 * Created on May 7, 2023, 9:20 PM
 */

#ifndef URWELLTOOLS_H
#define URWELLTOOLS_H
#include <vector>
#include <cmath>
#include <TH2D.h>
#include <TF1.h>

namespace uRwellTools {

    struct uRwellHit {
        int sector;
        int layer;
        int strip;
        int stripLocal;
        double adc;
        double adcRel;
        int ts;
        int slot;
    };

    
    const int clStripGap = 2; // The Max length of the gap in between strips with a given cluster 
    int getSlot(int ch);
    const int nSlot = 12; // The test Prototype has only 12 slots
    extern int slot_Offset[nSlot]; // Gives the 1st strip channel (unique channel) for the given slot
    
    const double strip_alpha = 10. * 0.017453293; // The strip angle in radians
    const double Y_0 = 250 + 723 * tan(strip_alpha);
    const double pitch = 1.; // mm
    
    const double uRWell_Y_max = 250.; // mm
    const double uRWell_Y_min = -250.; // mm
    const double uRwell_XTop = 728.;
    const double uRwell_XBot = 510.;

    /*
     * This function takes the h_in histogram as a first argument, following arguments are addresses of different 
     * efficiencies. It will calculates these efficiencies and write to corresponding addresses.
     * the "h_in" histogram is just the "number of V bv number of U" clusters.
     * eff_U is the efficiency of U cluster, eff_V is the efficiency of the V cluster, eff_OR is efficiency of "Any cluster", 
     * and eff_AND is the efficiency of having at least 1U and 1 V cluster
     */
    void CalcEfficiencies( TH2* h_in, double &eff_U, double &eff_V, double &eff_OR, double &eff_AND );
    
    /**
     * This function takes the h_in histogram as an input. It assumes the h_in histogram
     * is the cross "Y vs X" histogram. It projects it on "X" axis then fits with a Gaus + Pol4 function,
     * and calculates the number corresponding to the Gaussian.
     */
    double getNofBgrSbtrCrosses(TH2* h_in);
    
    /*
     * This function take cross x and y coordiantes, and checks whether the cross is inside the detector.
     * This assumes idealized geometry, and checks are based on variables uRWell_Y_max, uRWell_Y_min, uRwell_XTop, uRwell_XBot
     */
    bool IsInsideDetector( double x, double y );
    
    class uRwellCluster {
    public:
        uRwellCluster();

        void setEnergy(double aEnergy) {
            fEnergy = aEnergy;
        }
        
        void setHits(std::vector<uRwellHit>);
        std::vector<uRwellHit>* getHits(){
            return &fv_Hits;
        }
        void setNStrips(int aNStrips){
            fnStrips = fv_Hits.size();
        }

        const double getPeakADC(){
            return fPeakADC;
        }
        
        const double getAvgStrip(){
            return fAvgStrip;
        }
        
        // This method will be called when all hits of the cluster are identifed, and set into the cluster
        // I till find the Peak energy, the weighted avg strip etc.
        void FinalizeCluster();
        
    private:
        int fnStrips;
        double fEnergy;
        double fPeakADC;
        double fAvgStrip;
        std::vector<uRwellHit> fv_Hits;

        void findPeakEnergy();
        void findAvgStrip();

    };

    std::vector<uRwellCluster> getGlusters(std::vector<uRwellHit>);
    double getCrossX(double strip_U, double strip_V); // returns Cross UxV cross X coordinate
    double getCrossY(double strip_U, double strip_V); // returns Cross UxV cross Y coordinate
}

#endif /* URWELLTOOLS_H */