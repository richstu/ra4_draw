ra4_draw
========

Repository for plotting and table-making utilities combining the flexibility of analysis_code's simple and flexible string-based plotting with the speed of ra4_macros' "looping" scripts.

#### Commmands for RA4 analysis plots

To get just the pure MC kappa plots:
    
    ./compile.py && ./run/ra4/kappa_plots_tables.exe --mm off -m signal --mc_kappas -y 0 -n

Check kappa values in 5-6 jet CR, 2l+veto CR and SR+lowMET using MC in place of data:

    ./compile.py && ./run/ra4/kappa_plots_tables.exe --mm data -m cr56j --data_kappas -y 0
    ./compile.py && ./run/ra4/kappa_plots_tables.exe --mm data -m cr2lveto --data_kappas -y 0
    ./compile.py && ./run/ra4/kappa_plots_tables.exe --mm data -m lowmet --data_kappas -y 0

Unblinding SRs:
    
    ./compile.py && ./run/ra4/kappa_plots_tables.exe --mm data -m signal --data_kappas -y 0 -u

Additional options:

    * Use flag `--data_kappas` to make the kappa plots comparing mc to data or mc to mc w/ data stat
    * Use flag `--mc_kappas` to make the kappa plots showing just nominal MC kappas
    * Add the flag `--debug` (or `-d`) to see the cuts for each region. 
    * Use `--mm data` to compare actual data to MC in the CRs. Use `--mm off` to use MC as pseudodata.
    * Add the flag `--unblind` (or `-u`) to plot data in SRs as well.

For partial kappa plots in the systematic section:

    ./compile.py && ./run/ra4/plot_partial_kappa.exe -t lowmj -y 0 # plot all bkg low-MJ
    ./compile.py && ./run/ra4/plot_partial_kappa.exe -t highmj -y 0 # plot all bkg high-MJ
    ./compile.py && ./run/ra4/plot_partial_kappa.exe -t highmj -y 0 -n 1 -d 1 # plot wjets high-MJ
    ./compile.py && ./run/ra4/plot_partial_kappa.exe -t highmj -y 0 -n 2 -d 2 # plot ttjets high-MJ

Plots showing that CRs mimic behavior of SRs:

    ./compile.py && ./run/ra4/kappa_plots_tables.exe --mm mismeas_kappa -m cr56j_metbins --data_kappas -y 0
    ./compile.py && ./run/ra4/kappa_plots_tables.exe --mm mismeas_kappa -m cr2lveto_metbins --data_kappas -y 0
    ./compile.py && ./run/ra4/kappa_plots_tables.exe --mm mismeas_kappa -m signal_metbins --data_kappas -y 0
    ./compile.py && ./run/ra4/kappa_plots_tables.exe --mm w_isr -m cr2lveto_njbins --data_kappas -y 0
    ./compile.py && ./run/ra4/kappa_plots_tables.exe --mm w_isr -m signal_njbins --data_kappas -y 0

Generate a specific datacard, or wihtout signal uncertainties for the results plot:

    ./compile.py && ./run/ra4/write_datacards.exe -p 2100_100 -y 0 -u -d
    ./compile.py && ./run/ra4/write_datacards.exe -p 2100_100 -y 0 -u -d --no_syst

Getting the detailed signal systematics tables:

    ./python/sig_sys_table.py systable_lowmj.tex sys_SMS-T1tttt_mGluino-2100_mLSP-100_0_nom.txt sys_SMS-T1tttt_mGluino-1900_mLSP-1250_0_nom.txt --lowmj
    ./python/sig_sys_table.py systable_highmj.tex sys_SMS-T1tttt_mGluino-2100_mLSP-100_0_nom.txt sys_SMS-T1tttt_mGluino-1900_mLSP-1250_0_nom.txt

Get limits, this also generates the datacards for the full scan:
    
    ./python/send_limits.py

To perform the fit without the signal region observations, start by creating a workspace in order to create the masking variables:

    text2workspace.py datacard_SMS-T1tttt_mGluino-2100_mLSP-100_0_nom.txt --channel-masks

The resulting workspace can be examined using [dump_workspace.cxx](src/ra4/dump_workspace.cxx). Then, to do the fit with masked R4s:

    combine datacard_SMS-T1tttt_mGluino-2100_mLSP-100_0_nom.root -M FitDiagnostics --forceRecreateNLL --saveWorkspace --saveWithUncertainties --saveOverall --setParameters=mask_r4_lmet_lnb_lnj_lmj=1,mask_r4_lmet_lnb_lnj_hmj=1,mask_r4_lmet_lnb_hnj_lmj=1,mask_r4_lmet_lnb_hnj_hmj=1,mask_r4_lmet_mnb_lnj_lmj=1,mask_r4_lmet_mnb_lnj_hmj=1,mask_r4_lmet_mnb_hnj_lmj=1,mask_r4_lmet_mnb_hnj_hmj=1,mask_r4_lmet_hnb_lnj_lmj=1,mask_r4_lmet_hnb_lnj_hmj=1,mask_r4_lmet_hnb_hnj_lmj=1,mask_r4_lmet_hnb_hnj_hmj=1,mask_r4_mmet_lnb_lnj_lmj=1,mask_r4_mmet_lnb_lnj_hmj=1,mask_r4_mmet_lnb_hnj_lmj=1,mask_r4_mmet_lnb_hnj_hmj=1,mask_r4_mmet_mnb_lnj_lmj=1,mask_r4_mmet_mnb_lnj_hmj=1,mask_r4_mmet_mnb_hnj_lmj=1,mask_r4_mmet_mnb_hnj_hmj=1,mask_r4_mmet_hnb_lnj_lmj=1,mask_r4_mmet_hnb_lnj_hmj=1,mask_r4_mmet_hnb_hnj_lmj=1,mask_r4_mmet_hnb_hnj_hmj=1,mask_r4_hmet_lnb_lnj_lmj=1,mask_r4_hmet_lnb_lnj_hmj=1,mask_r4_hmet_lnb_hnj_lmj=1,mask_r4_hmet_lnb_hnj_hmj=1,mask_r4_hmet_mnb_lnj_lmj=1,mask_r4_hmet_mnb_lnj_hmj=1,mask_r4_hmet_mnb_hnj_lmj=1,mask_r4_hmet_mnb_hnj_hmj=1,mask_r4_hmet_hnb_lnj_lmj=1,mask_r4_hmet_hnb_lnj_hmj=1,mask_r4_hmet_hnb_hnj_lmj=1,mask_r4_hmet_hnb_hnj_hmj=1 --name=_nor4

Repeat w/o the `--setParameters` argument to get the global fit including R4's. The result plot can then be made using the output files from the above fit with the [python/plot_results.py](python/plot_results.py) script.

Goodness of fit estimate:

    combine -M GoodnessOfFit datacard_SMS-T1tttt_mGluino-2100_mLSP-100_0_nom.txt --algo=saturated --fixedSignalStrength=0
    combine -M GoodnessOfFit datacard_SMS-T1tttt_mGluino-2100_mLSP-100_0_nom.txt --algo=saturated --fixedSignalStrength=0 --toysFreq -t 100

Then check how many of the toys are above the value gotten from the data in the first command in the "limit" branch.


#### Code documentation
Doxygen-based documentation is available at [Adam's UCSB-HEP webpage](http://hep.ucsb.edu/people/ald77/documentation/doc_ra4_draw/).

#### Setup and compilation
Compilation requires c++11 and ROOT, but not CMSSW. To compile, simply run

    ./compile.py

#### Making histograms
An example script is available under src/test.cxx. To execute, compile and then run

    ./run/core/test.exe
