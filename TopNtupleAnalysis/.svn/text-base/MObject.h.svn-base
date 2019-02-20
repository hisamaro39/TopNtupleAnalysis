/**
 * @brief Base class for all Analysis objects (electrons, muons, jets, large-R jets, etc).
 * @author Danilo Enoque Ferreira de Lima <dferreir@cern.ch>
 */
#ifndef MOBJECT_H
#define MOBJECT_H

#include "TLorentzVector.h"
#include <string>
#include <map>

class MObject {

  public:

    enum Type {
      e = 11,
      mu = 13,
      jet = 1,
      part = 0,
      largejet = 99,
    };

    MObject();
    MObject(const TLorentzVector &v, const MObject::Type t = MObject::jet);
    virtual ~MObject();

    MObject::Type &type();
    const MObject::Type &type() const;
    TLorentzVector &mom();
    const TLorentzVector &mom() const;

    template <class C>
    float minDeltaR(const std::vector<C> &o, int &jidx) const {
      float dR = 99;
      jidx = -1;
      for (int k = 0; k < o.size(); ++k) {
        float tdR = mom().DeltaR(o[k].mom());
        if (tdR < dR) {
          dR = tdR;
          jidx = k;
        }
      }
      return dR;
    }

    template <class C>
    float minDeltaR(const std::vector<C> &o) const {
      float dR = 99;
      for (int k = 0; k < o.size(); ++k) {
        float tdR = mom().DeltaR(o[k].mom());
        if (tdR < dR) {
          dR = tdR;
        }
      }
      return dR;
    }

  protected:
    TLorentzVector m_mom;
    MObject::Type m_type;
    std::map<std::string, TLorentzVector> m_corrMap;
};

#endif

