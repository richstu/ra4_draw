ra4_draw
========

Repository for plotting and table-making utilities combining the flexibility of analysis_code's simple and flexible string-based plotting with the speed of ra4_macros' "looping" scripts.

#### Commmands for RA4 analysis plots

Check kappa values in 5-6 jet CR, 2l+veto CR and SR+lowMET using MC in place of data:

    ./compile.py && ./run/ra4/kappa_plots_tables.exe --mm off -m cr56j --data_kappas -y 0
    ./compile.py && ./run/ra4/kappa_plots_tables.exe --mm off -m cr2lveto --data_kappas -y 0
    ./compile.py && ./run/ra4/kappa_plots_tables.exe --mm off -m signal_met100 --mc_kappas -y 0 -d

Additional options:
    * Use flag `--data_kappas` to make the kappa plots comparing mc to data or mc to mc w/ data stat
    * Use flag `--mc_kappas` to make the kappa plots showing just nominal MC kappas
    * Add the flag `--debug` (or `-d`) to see the cuts for each region. 
    * Use `--mm data` to compare actual data to MC in the CRs.
    * Add the flag `--unblind` (or `-u`) to plot data in SRs as well.

For partial kappa plots in the systematic section:
    ./compile.py && ./run/ra4/plot_partial_kappa.exe -t lowmj -y 0 # plot all bkg low-MJ
    ./compile.py && ./run/ra4/plot_partial_kappa.exe -t highmj -y 0 # plot all bkg high-MJ
    ./compile.py && ./run/ra4/plot_partial_kappa.exe -t highmj -y 0 -n 1 -d 1 # plot wjets high-MJ
    ./compile.py && ./run/ra4/plot_partial_kappa.exe -t highmj -y 0 -n 2 -d 2 # plot ttjets high-MJ

#### Code documentation
Doxygen-based documentation is available at [Adam's UCSB-HEP webpage](http://hep.ucsb.edu/people/ald77/documentation/doc_ra4_draw/).

#### Setup and compilation
Compilation requires c++11 and ROOT, but not CMSSW. To compile, simply run

    ./compile.py

#### Making histograms
An example script is available under src/test.cxx. To execute, compile and then run

    ./run/core/test.exe
