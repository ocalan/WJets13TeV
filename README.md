# WJets13TeV
TreeAnalysis

Clone this file to either TreeAnalysis2012 or TreeAnalysis2011

For 2012 data analysis do: git clone git@github.com:iihe-cms-sw/TreeAnalysis.git TreeAnalysis2012

For 2011 data analysis do: git clone git@github.com:iihe-cms-sw/TreeAnalysis.git TreeAnalysis2011

if you want to run on: 1) WJets: move DYJets directory to WJets 2) TTbar to emu : move DYJets directory to TTbarEMu

Be aware that a rootlogon.C is provided. This is needed to load ROoUnfold automatically.

Once you have cloned the repository, you just have to build the RooUnfold-1.1.1 that is provided. This RooUnfold is a sloghtly modified version of the one you can download at CERN. That is why it is provided here.

To build it, nothing hard:

$ cd DYJets/RooUnfold-1.1.1 $ make

Bravo ! You have compiled the RooUnfold classes and you are ready to use it.

You can also create some directories that will be needed later:

$ mkdir HistoFiles PNGFiles

These first directory will contain the .root files created by runDYJets and co. and containing the histograms used during the whole process up to cross section measurements. The second will contain the .png and .pdf plots for all kind of comparison at several different levels.

These directories are however not considered by git (as can be seen in the .gitignore file).

Next you need to have the data to run on. For this, if your running on m-machine, a simple symbolic link will be enough. You can also opy the full directory.

Go to the racine directory TreeAnalysis2012 $ ln -s /user/tseva/analysis/TreeAnalysis2012/Data_Z_5311 Data_Z_5311

$ ln -s /user/aleonard/../tseva/analysis/TreeAnalysis2012/DataTTbarEMu DataTTbarEMu

$ ln -s /user/aleonard/../tseva/analysis/TreeAnalysis2012/DataW DataW
########## quick fix for lxplus from Darin
Set up CMSSW ( for root )

cmsrel CMSSW_5_3_0 cd CMSSW_5_3_0/src cmsenv
Get Git Repository

git clone git@github.com:iihe-cms-sw/TreeAnalysis.git TreeAnalysis2012
Setup W+Jets directory

cd TreeAnalysis2012 mkdir WJets mkdir WJets/DataW cp -r DYJets/* WJets/ cd WJets

--- Get example file cmsStage /store/group/phys_smp/WPlusJets/NtuplesTomislav/SMu_8TeV_T_s_channel_dR_5311.root DataW/SMu_8TeV_T_s_channel_dR_5311.root

--- Output directories mkdir HistoFiles PNGFiles

---- Compile RooUnfold cd RooUnfold-1.1.1 make cd -

--- Quick hacks for making LHAPDF work sed -i 's#/user/aleonard/LHAPDF/lib/libLHAPDF.so#/afs/cern.ch/cms/slc5_amd64_gcc434/external/lhapdf/5.8.5/lib/libLHAPDF.so#g' runDYJets.cc sed -i 's#/user/aleonard/lhapdf-5.9.1/include/#/afs/cern.ch/cms/slc5_amd64_gcc434/external/lhapdf/5.8.5/include/#g' rootlogon.C sed -i 's#NNPDF23_nlo_as_0118#NNPDF20_as_0118_100#g' ComputePDFUncertainties.cc sed -i 's#NNPDF23_nlo_as_0118#NNPDF20_as_0118_100#g' runDYJets.cc

--- Run the root file maker root -b runDYJets.cc
