#include "plot_maker.hpp"

#include "TCanvas.h"

#include "utilities.hpp"
#include "timer.hpp"
#include "thread_pool.hpp"
#include "named_func.hpp"

using namespace std;
using ScalarType = NamedFunc::ScalarType;
using VectorType = NamedFunc::VectorType;
using ScalarFunc = NamedFunc::ScalarFunc;
using VectorFunc = NamedFunc::VectorFunc;

namespace{
  bool HavePass(const VectorType &v){
    for(const auto&x: v){
      if(x) return true;
    }
    return false;
  }

  bool HavePass(const VectorType &a,
                const VectorType &b){
    for(auto it_a = a.cbegin(), it_b = b.cbegin();
        it_a != a.cend() && it_b != b.cend();
        ++it_a, ++it_b){
      if((*it_a) && (*it_b)) return true;
    }
    return false;
  }

  bool HavePass(const VectorType &a,
                const VectorType &b,
                const VectorType &c){
    for(auto it_a = a.cbegin(), it_b = b.cbegin(), it_c = c.cbegin();
        it_a != a.cend() && it_b != b.cend() && it_c != c.cend();
        ++it_a, ++it_b, ++it_c){
      if((*it_a) && (*it_b) && (*it_c)) return true;
    }
    return false;
  }
}

PlotMaker::PlotMaker():
  file_extensions_{"pdf"},
  stacks_(){
}

void PlotMaker::AddPlot(const vector<shared_ptr<Process> > &processes,
                        const HistoDef &histo_def){
  //Adds a plot to the list of plots to be produced. Does NOT fill or draw the histogram.
  stacks_.push_back(HistoStack(processes, histo_def));
}

void PlotMaker::MakePlots(){
  //Processes this list of plots provided with AddPlot and writes the results to disk
  if(file_extensions_.size() == 0){
    DBG("No file extensions provided => no plots produced.");
    return;
  }
  FillHistograms();
  PrintPlots();
}

void PlotMaker::FillHistograms(){
  //Iterates over all processes needed for requested plots and fills the necessary histograms
  set<shared_ptr<Process> > processes = GetProcesses();
  ThreadPool tp(thread::hardware_concurrency());
  for(const auto &proc: processes){
    tp.Push(bind(&PlotMaker::FillHistogram, this, ref(proc)));
  }
}

void PlotMaker::FillHistogram(const shared_ptr<Process> &proc){
  cout << "Filling histograms for the " << proc->name_ << " process..." << endl;

  Baby &baby = *(proc->baby_);
  
  vector<pair<HistoDef, TH1D * const> > histos;
  histos = GetHistos(proc);

  bool values_b, hist_cuts_b, proc_cuts_b, weights_b;
  ScalarType values_s, hist_cuts_s, proc_cuts_s, weights_s;
  VectorType values_v, hist_cuts_v, proc_cuts_v, weights_v;

  long num_entries = baby.GetEntries();
  Timer timer(num_entries, 1.);
  timer.Start();
  for(long entry = 0; entry < num_entries; ++entry){
    timer.Iterate();
    baby.GetEntry(entry);

    proc_cuts_b = proc->cut_.IsScalar();
    if(proc_cuts_b){
      proc_cuts_s = proc->cut_.GetScalar(baby);
      if(!proc_cuts_s) continue;
    }else{
      proc_cuts_v = proc->cut_.GetVector(baby);
    }
    if(!(proc_cuts_b || HavePass(proc_cuts_v))) continue;
      
    for(auto &histo: histos){
      const HistoDef &histo_def = histo.first;
      TH1D * const hist = histo.second;

      values_b = histo_def.var_.IsScalar();
      hist_cuts_b = histo_def.cut_.IsScalar();
      weights_b = histo_def.weight_.IsScalar();

      size_t min_vec_size = proc_cuts_b ? static_cast<size_t>(-1) : proc_cuts_v.size();
      if(hist_cuts_b){
        hist_cuts_s = histo_def.cut_.GetScalar(baby);
        if(!hist_cuts_s) continue;
      }else{
        hist_cuts_v = histo_def.cut_.GetVector(baby);
        if(hist_cuts_v.size() < min_vec_size) min_vec_size = hist_cuts_v.size();
      }
      if(!(proc_cuts_b || hist_cuts_b || HavePass(proc_cuts_v, hist_cuts_v))) continue;

      if(weights_b){
        weights_s = histo_def.weight_.GetScalar(baby);
        if(weights_s == 0.) continue;
      }else{
        weights_v = histo_def.weight_.GetVector(baby);
        if(weights_v.size() < min_vec_size) min_vec_size = weights_v.size();
      }
      if(!(proc_cuts_b || hist_cuts_b || weights_b || HavePass(proc_cuts_v, hist_cuts_v, weights_v))) continue;
      
      if(values_b){
        values_s = histo_def.var_.GetScalar(baby);
      }else{
        values_v = histo_def.var_.GetVector(baby);
        if(values_v.size() < min_vec_size) min_vec_size = values_v.size();
      }
      
      for(size_t i = 0;
          i < min_vec_size
            && (!(proc_cuts_b && hist_cuts_b && weights_b && values_b) || i < 1);
          ++i){
        ScalarType proc_cut = proc_cuts_b ? proc_cuts_s : proc_cuts_v.at(i);
        ScalarType hist_cut = hist_cuts_b ? hist_cuts_s : hist_cuts_v.at(i);
        if(!(proc_cut && hist_cut)) continue;
        ScalarType value = values_b ? values_s : values_v.at(i);
        ScalarType weight = weights_b ? weights_s : weights_v.at(i);
        hist->Fill(value, weight);
      }
    }
  }
}

void PlotMaker::PrintPlots(){
  //Takes already filled histograms and prints to file
  if(file_extensions_.size() == 0){
    DBG("No file extensions provided => no plots produced.");
    return;
  }
  TCanvas c;
  for(auto &stack: stacks_){
    bool drawn = false;
    vector<TH1D> bkgs(stack.backgrounds_.size());
    for(size_t ihist = stack.backgrounds_.size() - 1; ihist < stack.backgrounds_.size(); --ihist){
      bkgs.at(ihist) = stack.backgrounds_.at(ihist).second;
      for(size_t jhist = 0; jhist < ihist; ++jhist){
        bkgs.at(ihist) = bkgs.at(ihist) + stack.backgrounds_.at(jhist).second;
      }
      bkgs.at(ihist).Draw(drawn ? "hist same" : "hist");
      drawn = true;
    }
    for(auto &signal: stack.signals_){
      signal.second.Draw(drawn ? "hist same" : "hist");
      drawn = true;
    }
    for(auto &data: stack.datas_){
      data.second.Draw(drawn ? "hist same" : "hist");
      drawn = true;
    }

    string base_name = "plots/"+stack.definition_.GetName();
    for(const auto &ext: file_extensions_){
      string full_name = base_name+"."+ext;
      c.Print(full_name.c_str());
      cout << "Wrote plot to " << full_name << "." << endl;
    }
  }
}

set<shared_ptr<Process> > PlotMaker::GetProcesses() const{
  //Finds list of all processes needed for requested plots
  set<shared_ptr<Process> > processes;
  for(const auto &stack: stacks_){
    for(const auto &proc: stack.GetProcesses()){
      processes.insert(proc);
    }
  }
  return processes;
}

vector<pair<HistoDef, TH1D * const> > PlotMaker::GetHistos(const shared_ptr<Process> &process){
  //Gets list of plots that include a particular process
  vector<pair<HistoDef, TH1D * const> > histos;
  for(auto &stack: stacks_){
    auto procs = stack.GetProcesses();
    auto loc = procs.find(process);
    if(loc == procs.end()) continue;
    histos.push_back(pair<HistoDef, TH1D * const>(stack.definition_, stack.GetHisto(process)));
  }
  return histos;
}

vector<pair<HistoDef, const TH1D * const> > PlotMaker::GetHistos(const shared_ptr<Process> &process) const{
  //Gets list of plots that include a particular process
  vector<pair<HistoDef, const TH1D * const> > histos;
  for(const auto &stack: stacks_){
    auto procs = stack.GetProcesses();
    auto loc = procs.find(process);
    if(loc == procs.end()) continue;
    histos.push_back(pair<HistoDef, const TH1D * const>(stack.definition_, stack.GetHisto(process)));
  }
  return histos;
}
