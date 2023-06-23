#include <uRwellTools.h>
#include <algorithm>
#include <iostream>

using namespace std;

int uRwellTools::getSlot(int ch) {

    if (ch <= 64) {
        return 0;
    } else if (ch >= 65 && ch <= 192) {
        return 1;
    } else if (ch >= 193 && ch <= 320) {
        return 2;
    } else if (ch >= 321 && ch <= 448) {
        return 4;
    } else if (ch >= 449 && ch <= 576) {
        return 9;
    } else if (ch >= 577 && ch <= 704) {
        return 11;
    } else if (ch >= 1001 && ch <= 1064) {
        return 6;
    } else if (ch >= 1065 && ch <= 1192) {
        return 7;
    } else if (ch >= 1193 && ch <= 1320) {
        return 8;
    } else if (ch >= 1321 && ch <= 1448) {
        return 10;
    } else if (ch >= 1449 && ch <= 1576) {
        return 3;
    } else if (ch >= 1577 && ch <= 11704) {
        return 5;
    } else {
        return -1; // Should not happen, non existing slot
    }

}

int uRwellTools::slot_Offset[uRwellTools::nSlot] = {0, 64, 192, 1448, 320, 1576, 1000, 1064, 1192, 448, 1320, 576};

namespace uRwellTools {

    uRwellCluster::uRwellCluster() {
        fnStrips = 0;
        fEnergy = 0;
    }

    void uRwellCluster::setHits(std::vector<uRwellHit> aHits) {
        fv_Hits = aHits;
    }

    void uRwellCluster::findPeakEnergy() {

        double MaxADC = 0;

        for (auto curHit : fv_Hits) {
            MaxADC = MaxADC > curHit.adc ? MaxADC : curHit.adc;
        }

        fPeakADC = MaxADC;
    }

    void uRwellCluster::findAvgStrip() {

        double WeightedSum = 0;
        double ADCSum = 0;

        for (auto curHit : fv_Hits) {
            WeightedSum = WeightedSum + curHit.adc * curHit.strip;
            ADCSum = ADCSum + curHit.adc;
        }

        fAvgStrip = WeightedSum / ADCSum;
    }

    void uRwellCluster::FinalizeCluster() {
        findPeakEnergy();
        findAvgStrip();
    }

    std::vector<uRwellCluster> getGlusters(std::vector<uRwellHit> v_Hits) {
        //
        // Let's 1st sort the vector of strips
        //
        vector<uRwellCluster> v_Clusters;
        sort(v_Hits.begin(), v_Hits.end(), [ ](const auto& lhs, const auto& rhs) {
            return lhs.strip < rhs.strip;
        });


        double clEnergy = 0;
        bool newCluster = true;
        int prev_Strip = -10000; // Some number that clearly is not a real strip number
        //        vector<int> v_strips;
        vector<uRwellHit> v_ClHits;

        for (int i = 0; i < v_Hits.size(); i++) {

            int curStrip = v_Hits.at(i).strip;
            double curHitEnergy = v_Hits.at(i).adc;

            //cout << "v_Hits.at(i).strip = " << v_Hits.at(i).strip << "    curStrip =  " << curStrip << endl;
            if ( (curStrip - prev_Strip <= (clStripGap + 1)) || i == 0) {
                clEnergy = clEnergy + curHitEnergy;
                v_ClHits.push_back(v_Hits.at(i));

                if (i == v_Hits.size() - 1) {
                    uRwellTools::uRwellCluster curCluster;
                    curCluster.setEnergy(clEnergy);
                    curCluster.setHits(v_ClHits);
                    curCluster.FinalizeCluster();
                    v_Clusters.push_back(curCluster);
                }

            } else {
                uRwellTools::uRwellCluster curCluster;
                curCluster.setEnergy(clEnergy);
                curCluster.setHits(v_ClHits);
                curCluster.FinalizeCluster();
                v_Clusters.push_back(curCluster);

                v_ClHits.clear();
                v_ClHits.shrink_to_fit();

                clEnergy = curHitEnergy;
                v_ClHits.push_back(v_Hits.at(i));
            }
            prev_Strip = curStrip;
        }

        //cout<<"The size of the cluster is "<<v_Clusters.size()<<endl;
        return v_Clusters;
    }

    double getCrossX(double strip_U, double strip_V) {
        return pitch * (strip_U - strip_V) / (2 * sin(strip_alpha));
    }

    double getCrossY(double strip_U, double strip_V) {
        double crs_x = getCrossX(strip_U, strip_V);
        return tan(strip_alpha) * crs_x + Y_0 - (strip_U * pitch) / cos(strip_alpha);
    }

}