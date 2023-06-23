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