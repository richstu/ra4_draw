#include <iostream>
#include <iomanip>

#include "TFile.h"
#include "TIterator.h"

#include "RooArgSet.h"
#include "RooRealVar.h"

#include "RooStats/ModelConfig.h"

using namespace std;
using namespace RooStats;

int main(int argc, char *argv[]){
  for(int iarg = 1; iarg < argc; ++iarg){
    TFile file(argv[iarg], "read");
    if(!file.IsOpen()) continue;
    RooWorkspace *w = static_cast<RooWorkspace*>(file.Get("w"));
    if(w == nullptr) continue;
    w->Print();

    cout << fixed << setprecision(2);

    cout << "Constants:" << endl;
    const RooArgSet &vars = w->allVars();
    TIterator *iter_ptr = vars.createIterator();
    for(; iter_ptr != nullptr && *(*iter_ptr) != nullptr; iter_ptr->Next()){
      RooAbsArg *var = static_cast<RooAbsArg*>(*(*iter_ptr));
      if(var == nullptr) continue;
      if(!var->isConstant()) continue;
      var->Print();
    }

    cout << endl << "Fundamental variables:" << endl;
    iter_ptr->Reset();
    for(; iter_ptr != nullptr && *(*iter_ptr) != nullptr; iter_ptr->Next()){
      RooAbsArg *var = static_cast<RooAbsArg*>(*(*iter_ptr));
      if(var == nullptr) continue;
      if(var->isConstant()) continue;
      var->Print();
    }

    cout << endl << "Yields: " << endl;
    cout
      << ' ' << setw(64) << "Name"
      << ' ' << setw(12) << "Value"
      << endl;

    iter_ptr->Reset();
    for(; iter_ptr != nullptr && *(*iter_ptr) != nullptr; iter_ptr->Next()){
      RooAbsArg *var = static_cast<RooAbsArg*>(*(*iter_ptr));
      if(var == nullptr) continue;
      string nobs_name = var->GetName();
      RooRealVar *nobs_arg = static_cast<RooRealVar*>(w->arg(nobs_name.c_str()));
      double nobs = nobs_arg->getVal();
      cout
    << ' ' << setw(64) << nobs_name
    << ' ' << setw(12) << nobs
    << endl;
    }
    
  }
}
