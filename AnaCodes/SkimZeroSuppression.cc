/* 
 * File:   SkimZeroSuppression.cc
 * Author: rafopar
 *
 * Created on May 4, 2023, 4:03 PM
 */

#include <cstdlib>

#include <TH2D.h>
#include <TH1D.h>
#include <TMath.h>
#include <TFile.h>

// ===== Hipo headers =====
#include <reader.h>
#include <writer.h>
#include <dictionary.h>

#include <uRwellTools.h>

using namespace std;
int getSlot(int uniqueChannel);

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

/*
 * 
 */
int main(int argc, char** argv) {

    char outputFile[256];
    char inputFile[256];

    int run = 0;
    int fnum = -1;
    if (argc > 2) {
        run = atoi(argv[1]);
        fnum = atoi(argv[2]);
        //sprintf(inputFile, "%s", argv[1]);
        sprintf(inputFile, "Data/decoded_%d_%d.hipo", run, fnum);
        //sprintf(inputFile, "Data/decoded_%d.hipo", run);
        sprintf(outputFile, "Skim_ZeroSuppr_%d_%d.hipo", run, fnum);
    } else {
        std::cout << " *** please provide a run number..." << std::endl;
        exit(0);
    }

    /**
     * Creating schema for uRwell::Hit
     */

    hipo::schema sch("uRwell::Hit", 90, 1);
    sch.parse("sec/S,layer/S,strip/S,stripLocal/S,adc/F,adcRel/F,ts/S,slot/S");
    sch.show();


    hipo::reader reader;
    reader.open(inputFile);

    hipo::dictionary factory;

    reader.readDictionary(factory);

    factory.show();
    hipo::event event;
    int evCounter = 0;

    hipo::bank buRWellADC(factory.getSchema("URWELL::adc"));
    hipo::bank bRAWADc(factory.getSchema("RAW::adc"));
    hipo::bank bRunConf(factory.getSchema("RUN::config"));

    hipo::writer writer;
    writer.getDictionary().addSchema(sch);
    writer.getDictionary().addSchema( factory.getSchema("RAW::adc") );
    writer.getDictionary().addSchema( factory.getSchema("RUN::config") );
    writer.open(outputFile);

    int __bank_Sec_INDEX_ = buRWellADC.getSchema().getEntryOrder("sector");
    int __bank_Layer_INDEX_ = buRWellADC.getSchema().getEntryOrder("layer");
    int __bank_Component_INDEX_ = buRWellADC.getSchema().getEntryOrder("component");
    int __bank_Order_INDEX_ = buRWellADC.getSchema().getEntryOrder("order");
    int __bank_ADC_INDEX_ = buRWellADC.getSchema().getEntryOrder("ADC");
    int __bank_Time_INDEX_ = buRWellADC.getSchema().getEntryOrder("time");
    int __bank_Ped_INDEX_ = buRWellADC.getSchema().getEntryOrder("ped");

    const double sigm_threshold = 3.; // Represents the threshold of the ADC in units of the sigma.
    const double sigm_thresholduRwell = 3.; // Represents the threshold of the ADC in units of the sigma.
    const double sigm_thresholduRwellMAX = 100.; // This is to cut noisy channels, that produce 


    const int nGEMChannels = 256;
    const int n_ts = 9;
    const int sec_GEM = 8; // GEM is in Sec 8
    const int sec_uRWell = 6; // uRWELL is in Sec 8


    /**
     * Reading the pedestal file and and fill a map for pedestals and RMSs
     */

    std::map<int, double> m_ped_mean; // Mean value of the pedestal
    std::map<int, double> m_ped_rms; // rms of the pedestal
    std::map<int, double> m_ped_GEM_mean; // Mean value of the pedestal
    std::map<int, double> m_ped_GEM_rms; // rms of the pedestal

    ifstream inp_ped(Form("PedFiles/Peds_%d", run));

    if (!inp_ped.is_open()) {
        cout << "Can not find the pedestal file \"PedFiles/CosmicPeds.dat\". " << endl;
        cout << "Exiting..." << endl;
        exit(1);
    }

    cout << "Kuku" << endl;
    while (!inp_ped.eof()) {
        int ch;
        double mean, rms;
        inp_ped>>ch;
        inp_ped>>mean;
        inp_ped>>rms;

        m_ped_mean[ch] = mean;
        m_ped_rms[ch] = rms;
    }

    ifstream inp_ped_GEM(Form("PedFiles/GEM_Peds_%d", run));

    if (!inp_ped_GEM.is_open()) {
        cout << "Can not find the pedestal file. " << endl;
        cout << "Exiting..." << endl;
        exit(1);
    }

    while (!inp_ped_GEM.eof()) {
        int ch;
        double mean, rms;
        inp_ped_GEM>>ch;
        inp_ped_GEM>>mean;
        inp_ped_GEM>>rms;

        m_ped_GEM_mean[ch] = mean;
        m_ped_GEM_rms[ch] = rms;
    }

    cout << "The pedestal map is loaded." << endl;

    try {

        while (reader.next() == true) {
            reader.read(event);

            evCounter = evCounter + 1;

            //if( evCounter > 200 ){break;}
            if (evCounter % 1000 == 0) {
                cout.flush() << "Processed " << evCounter << " events \r";
            }

            event.getStructure(buRWellADC);
            event.getStructure(bRAWADc);
            event.getStructure(bRunConf);

            int ev_Number = bRunConf.getInt("event", 0);

            int n_uRwellADC = buRWellADC.getRows();

            if (n_uRwellADC == 0) {
                continue;
            }


            double ADC_GEM_[nGEMChannels] = {0};
            double ADCRel_GEM_[nGEMChannels] = {0};

            std::map<int, double> m_ADC_uRWELL;
            std::map<int, double> m_ADCRel_uRWELL;


            for (int i = 0; i < n_uRwellADC; i++) {
                int sector = buRWellADC.getInt(__bank_Sec_INDEX_, i);
                int layer = buRWellADC.getInt(__bank_Layer_INDEX_, i);
                int channel = buRWellADC.getInt(__bank_Component_INDEX_, i);
                int ADC = buRWellADC.getInt(__bank_ADC_INDEX_, i);
                int uniqueChan = int(buRWellADC.getFloat(__bank_Time_INDEX_, i));
                int ts = buRWellADC.getInt(__bank_Ped_INDEX_, i);

                int slot = layer;

                if (sector == sec_GEM) {
                    //                    double ADC_OffsetCorrected = m_ped_GEM_mean[channel] - ADC;
                    //                    ADC_GEM_[channel-1] = ADC_GEM_[channel-1] + ADC_OffsetCorrected;
                    ADC_GEM_[channel - 1] = ADC_GEM_[channel - 1] + ADC;

                    //if (ADC == 0) {
                    //cout<<"Event# = "<<ev_Number<<"   channel = "<<channel<<"  ts = "<<ts << "   ADC = " << ADC << endl;
                    //}
                } else if (sector == sec_uRWell) {
                    m_ADC_uRWELL[uniqueChan] = m_ADC_uRWELL[uniqueChan] + double(ADC);
                }

            }


            vector<uRwellHit> v_GEM_Hits;
            for (int ich = 0; ich < nGEMChannels; ich++) {
                ADC_GEM_[ich] = m_ped_GEM_mean[ich + 1] - ADC_GEM_[ich] / double(n_ts);
                ADCRel_GEM_[ich] = ADC_GEM_[ich] / m_ped_GEM_rms[ich + 1];

                if (ADCRel_GEM_[ich] > sigm_threshold) {
                    uRwellHit curHit;
                    curHit.adc = ADC_GEM_[ich];
                    curHit.adcRel = ADCRel_GEM_[ich];
                    curHit.sector = sec_GEM;
                    curHit.layer = 0;
                    curHit.slot = 12 + ich / 128;
                    curHit.strip = ich;
                    curHit.stripLocal = ich%128; // For GEM this doesn't have a good meaning without knowing the strip to connector internal mapping
                    curHit.ts = 0; // For now we will not prserve the hit time information
                    v_GEM_Hits.push_back(curHit);
                }
            }


            vector<uRwellHit> v_uRwell_Hits;
            for (auto it = m_ADC_uRWELL.begin(); it != m_ADC_uRWELL.end(); ++it) {
                int ch = it->first;
                m_ADC_uRWELL[ch] = m_ped_mean[ch] - m_ADC_uRWELL[ch] / double(n_ts);
                m_ADCRel_uRWELL[ch] = m_ADC_uRWELL[ch] / m_ped_rms[ch];

                if (m_ADCRel_uRWELL[ch] > sigm_thresholduRwell) {
                    uRwellHit curHit;
                    curHit.adc = m_ADC_uRWELL[ch];
                    curHit.adcRel = m_ADCRel_uRWELL[ch];
                    curHit.sector = sec_uRWell;
                    curHit.layer = 1 + ch / 1000;
                    curHit.slot = getSlot(ch);
                    curHit.strip = ch % 1000;
                    curHit.stripLocal = ch - uRwellTools::slot_Offset[curHit.slot];
                    curHit.ts = 0;
                    v_uRwell_Hits.push_back(curHit);
                }
            }

            int n_TotHits = v_GEM_Hits.size() + v_uRwell_Hits.size();

            if (n_TotHits == 0) {
                continue;
            }

            hipo::bank buRwellHits(sch, n_TotHits);
            hipo::event outEvent;

            int col = 0;
            //     *********** Writing GEM hits above the threshold **********
            //sch.parse("sec/S,layer/S,strip/S,adc/F,adcRel/F,ts/S,slot/S");
            for (auto curHit : v_uRwell_Hits) {
                buRwellHits.putShort("sec", col, short(curHit.sector));
                buRwellHits.putShort("layer", col, short(curHit.layer));
                buRwellHits.putShort("strip", col, short(curHit.strip));
                buRwellHits.putShort("stripLocal", col, short(curHit.stripLocal));
                buRwellHits.putFloat("adc", col, float(curHit.adc));
                buRwellHits.putFloat("adcRel", col, float(curHit.adcRel));
                buRwellHits.putShort("ts", col, short(curHit.ts));
                buRwellHits.putShort("slot", col, short(curHit.slot));

                col = col + 1;
            }

            //     *********** Writing uRwell hits above the threshold **********
            for (auto curHit : v_GEM_Hits) {
                buRwellHits.putShort("sec", col, short(curHit.sector));
                buRwellHits.putShort("layer", col, short(curHit.layer));
                buRwellHits.putShort("strip", col, short(curHit.strip));
                buRwellHits.putShort("stripLocal", col, short(curHit.stripLocal));
                buRwellHits.putFloat("adc", col, float(curHit.adc));
                buRwellHits.putFloat("adcRel", col, float(curHit.adcRel));
                buRwellHits.putShort("ts", col, short(curHit.ts));
                buRwellHits.putShort("slot", col, short(curHit.slot));

                col = col + 1;
            }

            outEvent.addStructure(bRAWADc);
            outEvent.addStructure(bRunConf);
            outEvent.addStructure(buRwellHits);
            writer.addEvent(outEvent);

        }
    } catch (const char msg) {
        cerr << msg << endl;
    }

    writer.close();
    writer.showSummary();

    return 0;
}

int getSlot(int ch) {

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