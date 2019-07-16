#ifndef __AT_QC_HH__
#define __AT_QC_HH__

#include <vector>
#include "AT_Charged.h"

class TProfile;

class AT_QC : public AnalysisTask {
 public:
  AT_QC();
  virtual ~AT_QC();
  virtual void Init();
  virtual void Exec();
  virtual void Finish();

 private:
  TProfile *hd[4];
  TProfile *hr[4];
};

#endif
