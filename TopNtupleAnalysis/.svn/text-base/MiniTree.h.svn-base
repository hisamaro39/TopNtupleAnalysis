/**
 * @brief Class that reads information from the input file and puts it in a object-oriented format in Event.
 * @author Danilo Enoque Ferreira de Lima <dferreir@cern.ch>
 */
#ifndef MINITREE_H
#define MINITREE_H

#include <string>
#include "TChain.h"
#include "TTree.h"
#include <vector>
#include "TopNtupleAnalysis/Event.h"

#include "TFile.h"
#include <map>

using namespace std;

class MiniTree {
  public:

    MiniTree(bool toWrite = true, const std::string &file = "tree.root", const std::string &name="nominal");
    virtual ~MiniTree();

    void addFileToRead(const std::string &fname);
    void addFileToRead(const std::string &fname, const std::string &treeName);
    int GetEntries();
    double getSumWeights();
    void read(int event, Event &e);

    double &sumWeights();
    TTree *m_chain;
    double m_sumWeights;

    TFile  *m_file;
    TTree *m_num;
    std::string m_name;

    unsigned int         &ui(const std::string &n);
    ULong64_t            &ul64(const std::string &n);
    int                  &i(const std::string &n);
    float                &f(const std::string &n);
    char                 &c(const std::string &n);
    std::vector<char>    *vc(const std::string &n);
    std::vector<int>     *vi(const std::string &n);
    std::vector<float>   *vf(const std::string &n);
    std::vector<std::vector<int> > *vvi(const std::string &n);
    std::vector<std::vector<float> > *vvf(const std::string &n);

    enum MTType {mtUint = 0, mtInt, mtFloat, mtChar, mtVInt, mtVFloat, mtVChar, mtVVInt, mtVVFloat, mtULong64};
    std::map<std::string, MTType> m_brs;
    
  private:

    void prepareBranches();

    std::map<std::string, float> m_f;
    std::map<std::string, unsigned int> m_ui;
    std::map<std::string, ULong64_t> m_ul64;
    std::map<std::string, int> m_i;
    std::map<std::string, char> m_c;
    std::map<std::string, std::vector<char> *> m_vc;
    std::map<std::string, std::vector<float> *> m_vf;
    std::map<std::string, std::vector<int> *> m_vi;
    std::map<std::string, std::vector<std::vector<int> > *> m_vvi;
    std::map<std::string, std::vector<std::vector<float> > *> m_vvf;

};

#endif

