#ifndef SOLVER
#define SOLVER


const double PI=3.14159265358979323846264338;
const double TWOPI=2.*PI;
const double SQRTTHREEOVERTWO = 0.86602540378443864676;

// Use Eigen to diagonalize? Or LAPACK? Try Eigen.
#include <Eigen/Dense>
using namespace Eigen;

#include "Quantities.hpp"
#include "makepretty.h"

#include <fstream>


class Solver
{
public:
  string dir;
  int Nsites, TWOL, Nh, Ns, nu;
  int maxIndexValue;
  int maxStateValue;
  unsigned int twomax; //Are these used in  Solver?
  vector<int> TWOLpow; //Are these used in  Solver?

  double tl, tr, Jzl, Jzr, Jpml, Jpmr; //t, J_z and J_{\pm} for legs and rungs.

  indexstate converttable;

  Matrix<double,Dynamic,Dynamic> H;

  Eigen::Matrix<double, -1, 1, 0, -1, 1> eigenvals;
  Matrix<double,Dynamic,Dynamic> eigenvecs;

  bool OBC;

  bool PYROCHLORE;

  int CUTOFF;

  bool EIGVEC, ZEROCORR, NNCORR, MIDCORR;

  complex<double> zero;
  int Nb, Nt;

  ofstream Qout;


  Solver();
  //Solver(string dir0, int L0, int Nh0, double tl0, double tr0, double Jzl0, double Jzr0, double Jpml0, double Jpmr0, bool OBC0, bool PYROCHLORE0, bool EIGVEC0, bool CORR0);
  Solver(string dir0, ReadInputFiles params);
  void solve();

  void makebasis();

  unsigned int statevec_to_statenum(vector<short int> statevec);
  vector<short int> statenum_to_statevec(unsigned int statenum);

  void fillH();
  void diagonalise();

  Matrix<double,Dynamic,Dynamic> NhMATRIX(int siteind);
  vector<vector<complex<double>>> NhExp(vector<double> beta, vector<double> time);
  vector<vector<vector<complex<double>>>> NhCorrZERO(vector<double> beta, vector<double> time);
  vector<vector<vector<complex<double>>>> NhCorrNN(vector<double> beta, vector<double> time);
  vector<vector<vector<complex<double>>>> NhCorrMID(vector<double> beta, vector<double> time);

  vector<double> Szmat(int siteind);

  Matrix<double, Dynamic, Dynamic> SzmatMATRIX(int siteind);
  //complex<double> SzCorr(int i, int j, double beta, double t);

  vector<vector<vector<complex<double>>>> SzCorrMatZERO(vector<double> beta, vector<double> time);
  vector<vector<vector<complex<double>>>> SzCorrMatNN(vector<double> beta, vector<double> time);
  vector<vector<vector<complex<double>>>> SzCorrMatMID(vector<double> beta, vector<double> time);
  vector<vector<vector<complex<double>>>> SzCorr(vector<double> beta, vector<double> time);
  Matrix<double, Dynamic, Dynamic> makeSminus(indexstate converttablep);
  Matrix<double, Dynamic, Dynamic> SminusMATRIX(indexstate converttablep, int siteind);
  //complex<double> SpmCorr(int i, int j, Matrix<double, Dynamic, Dynamic> Sminusmat, Eigen::Matrix<double, -1, 1, 0, -1, 1> eigenvalsp, Matrix<double,Dynamic,Dynamic> eigenvecsp, indexstate converttablep, double beta, double t);
  vector<vector<vector<complex<double>>>> SpmCorrMatZERO(Eigen::Matrix<double, -1, 1, 0, -1, 1> eigenvalsp, Matrix<double,Dynamic,Dynamic> eigenvecsp, indexstate converttablep, vector<double> beta, vector<double> time);
  vector<vector<vector<complex<double>>>> SpmCorrMatNN(Eigen::Matrix<double, -1, 1, 0, -1, 1> eigenvalsp, Matrix<double,Dynamic,Dynamic> eigenvecsp, indexstate converttablep, vector<double> beta, vector<double> time);
  vector<vector<vector<complex<double>>>> SpmCorrMatMID(Eigen::Matrix<double, -1, 1, 0, -1, 1> eigenvalsp, Matrix<double,Dynamic,Dynamic> eigenvecsp, indexstate converttablep, vector<double> beta, vector<double> time);
  vector<vector<vector<complex<double>>>> SpmCorr(Matrix<double, Dynamic, Dynamic> Sminusmat, Eigen::Matrix<double, -1, 1, 0, -1, 1> eigenvalsp, Matrix<double,Dynamic,Dynamic> eigenvecsp, indexstate converttablep, vector<double> beta, vector<double> time);

  double Sx2(Eigen::Matrix<double, -1, 1, 0, -1, 1> statecoeffs);
  double Sy2(Eigen::Matrix<double, -1, 1, 0, -1, 1> statecoeffs);
  double Sz2(Eigen::Matrix<double, -1, 1, 0, -1, 1> statecoeffs);
  double S(Eigen::Matrix<double, -1, 1, 0, -1, 1> statecoeffs);

  vector<double> HoleDens(Eigen::Matrix<double, -1, 1, 0, -1, 1> statecoeffs);

  vector<short int> translate(vector<short int> statevec, int trans);
  vector<double> TransMat(int trans);
  void GSQvec();

  vector<short int> parity(vector<short int> statevec);
  Matrix<double, Dynamic, Dynamic> ParMat();
  void GSparity();

  string index_to_string(unsigned short int stateind);

  void PrintState(Eigen::Matrix<complex<double>, -1, 1, 0, -1, 1> state, ostream& os);

  void WriteEigvals();
  void WriteSzStot();
  void WritePartition(vector<double> beta, vector<double> partition);
  void WriteCorr(vector<double> beta, vector<double> time, vector<vector<vector<complex<double>>>> z, vector<vector<vector<complex<double>>>> pm, string file);
  void WriteHoleDens();
  void WriteHoleCorr(vector<double> beta, vector<double> time, vector<vector<vector<complex<double>>>> NhNh, string file);
  void resetdatafiles();

};

Solver::Solver(){}

Solver::Solver(string dir0, ReadInputFiles params)
{
  zero.real(0); zero.imag(0);

  Nb = 0; Nt =0;

  dir = dir0;

  Nsites = params.Nsites;
  TWOL = Nsites;
  Nh = params.Nh;
  nu = 0; //Default

  tl = params.tl;
  tr = params.tr;
  Jzl = params.Jzl;
  Jzr = params.Jzr;
  Jpml = params.Jpml;
  Jpmr = params.Jpmr;

  OBC = params.OBC;

  PYROCHLORE = params.PYROCHLORE;

  CUTOFF = params.CUTOFF;

  EIGVEC = params.EIGVECS;
  ZEROCORR = params.ZEROCORR;
  NNCORR = params.NNCORR;
  MIDCORR = params.MIDCORR;

  Ns = TWOL - Nh;
  twomax = 1<<Ns;
  TWOLpow = vector<int>(Nh);
  for (int i = 0; i < Nh; i++)
  {
    TWOLpow[i] = pow(TWOL, i);
  }

  if (params.RESETFILES)
  {
    Qout = ofstream(dir + "Qstates.txt");
    if (!Qout.is_open())
      cout<<"Could not open file" << endl;
  }
  else
  {
    Qout = ofstream(dir + "Qstates.txt", std::ios_base::app);
    if (!Qout.is_open())
      cout<<"Could not open file" << endl;
  }
}


void Solver::solve()
{
  vector<double> mineigvals(Ns+1);   //Minimum eigenvalue for each sector
  vector<double> mineigvalspm(Ns+1); //Minimum eigenvalues for a sector and the previous so that Corrpm doesn't blow up at large beta

  Eigen::Matrix<double, -1, 1, 0, -1, 1> eigenvalsp;
  Matrix<double,Dynamic,Dynamic> eigenvecsp;
  indexstate converttablep;

  vector<double> beta = {0, 0.5, 1, 2, 10, 50, 90, 99, 100, 200, 500, 1000, 10000};
  vector<double> time = {0};

  Nb = beta.size();
  Nt = time.size();

  double start, stop;

  vector<vector<double>> partfunc(Ns+1, vector<double>(Nb, 0.0));       //Partition function for each magnetisation sector, scaled by e^(-beta mineigvals)
  vector<vector<vector<vector<complex<double>>>>> corrfunczzZERO(Ns+1, vector<vector<vector<complex<double>>>>(TWOL, vector<vector<complex<double>>>(Nb, vector<complex<double>>(Nt, zero))));       //Correlation functions for each magnetisation sector.
  vector<vector<vector<vector<complex<double>>>>> corrfuncpmZERO(Ns+1, vector<vector<vector<complex<double>>>>(TWOL, vector<vector<complex<double>>>(Nb, vector<complex<double>>(Nt, zero))));       //Correlation functions for each magnetisation sector.
  vector<vector<vector<vector<complex<double>>>>> corrfuncNhZERO(Ns+1, vector<vector<vector<complex<double>>>>(TWOL, vector<vector<complex<double>>>(Nb, vector<complex<double>>(Nt, zero))));       //Correlation functions for each magnetisation sector.
  vector<vector<vector<vector<complex<double>>>>> corrfunczzNN(Ns+1, vector<vector<vector<complex<double>>>>(TWOL, vector<vector<complex<double>>>(Nb, vector<complex<double>>(Nt, zero))));       //Correlation functions for each magnetisation sector.
  vector<vector<vector<vector<complex<double>>>>> corrfuncpmNN(Ns+1, vector<vector<vector<complex<double>>>>(TWOL, vector<vector<complex<double>>>(Nb, vector<complex<double>>(Nt, zero))));       //Correlation functions for each magnetisation sector.
  vector<vector<vector<vector<complex<double>>>>> corrfuncNhNN(Ns+1, vector<vector<vector<complex<double>>>>(TWOL, vector<vector<complex<double>>>(Nb, vector<complex<double>>(Nt, zero))));       //Correlation functions for each magnetisation sector.

  vector<vector<vector<vector<complex<double>>>>> corrfunczzMID(Ns+1, vector<vector<vector<complex<double>>>>(TWOL, vector<vector<complex<double>>>(Nb, vector<complex<double>>(Nt, zero))));       //Correlation functions for each magnetisation sector.
  vector<vector<vector<vector<complex<double>>>>> corrfuncpmMID(Ns+1, vector<vector<vector<complex<double>>>>(TWOL, vector<vector<complex<double>>>(Nb, vector<complex<double>>(Nt, zero))));       //Correlation functions for each magnetisation sector.
  vector<vector<vector<vector<complex<double>>>>> corrfuncNhMID(Ns+1, vector<vector<vector<complex<double>>>>(TWOL, vector<vector<complex<double>>>(Nb, vector<complex<double>>(Nt, zero))));       //Correlation functions for each magnetisation sector.

  //Do nu=0 first:
  nu = 0;
  makebasis();
  fillH();
  //cout << H << endl;
  diagonalise();

  mineigvals[0] = eigenvals[0];
  mineigvalspm[0] = eigenvals[0];
  for(int b = 0; b < Nb; b++) partfunc[0][b] = partitionfunction(eigenvals, beta[b]);

  eigenvalsp = eigenvals;
  eigenvecsp = eigenvecs;
  converttablep = converttable;

  GSQvec();

  WriteEigvals();
  WriteSzStot();
  WriteHoleDens();
  if(ZEROCORR)
  {
    cout << "IN CORR FOR nu=0" << endl;
    cout << "BEFORE SZCORR FOR nu=0" << endl;

    corrfunczzZERO[0] = SzCorrMatZERO(beta, time); //Only contribution from Sz?
    for(int j = 0; j < TWOL; j++) for(int b = 0; b < Nb; b++) for(int t = 0; t < Nt; t++) corrfuncpmZERO[0][j][b][t] = zero;

    if (Nh > 1)
    {
      corrfuncNhZERO[0] = NhCorrZERO(beta, time);
    }
  }
  if(NNCORR)
  {
    cout << "IN CORR FOR nu=0" << endl;
    cout << "BEFORE SZCORR FOR nu=0" << endl;

    corrfunczzNN[0] = SzCorrMatNN(beta, time); //Only contribution from Sz?
    for(int j = 0; j < TWOL; j++) for(int b = 0; b < Nb; b++) for(int t = 0; t < Nt; t++) corrfuncpmNN[0][j][b][t] = zero;

    if (Nh > 1)
    {
      corrfuncNhNN[0] = NhCorrNN(beta, time);
    }
  }
  if(MIDCORR)
  {
    cout << "IN CORR FOR nu=0" << endl;
    cout << "BEFORE SZCORR FOR nu=0" << endl;

    corrfunczzMID[0] = SzCorrMatMID(beta, time); //Only contribution from Sz?
    for(int j = 0; j < TWOL; j++) for(int b = 0; b < Nb; b++) for(int t = 0; t < Nt; t++) corrfuncpmMID[0][j][b][t] = zero;

    if (Nh > 1)
    {
      corrfuncNhMID[0] = NhCorrMID(beta, time);
    }
  }

  Matrix<double, Dynamic, Dynamic> Sminus;

  for (int mynu = 1; mynu < Ns+1; mynu++)
  {
    nu = mynu;
    makebasis();
    fillH();
    //cout << H << endl;
    diagonalise();
    //cout << "nu = " << nu << endl;
    //cout << eigenvecs << endl;
    mineigvals[mynu] = eigenvals[0];
    mineigvalspm[mynu] = (eigenvalsp[0] < eigenvals[0]) ? eigenvalsp[0] : eigenvals[0];

    for(int b = 0; b < Nb; b++) partfunc[mynu][b] = partitionfunction(eigenvals, beta[b]);

    WriteEigvals();
    WriteSzStot();
    WriteHoleDens();

    GSQvec();

    //Compute correlations here!

    if(ZEROCORR)
    {
      Sminus = makeSminus(converttablep);

      cout << "BEFORE SZCORR FOR nu=" << nu << endl;

      start = clock(); //start clock
      corrfunczzZERO[mynu] = SzCorrMatZERO(beta, time);
      stop = clock();

      cout << "SzMAT time for nu = " << mynu << ": " << (stop-start)/CLOCKS_PER_SEC << endl;

      start = clock(); //start clock
      corrfuncpmZERO[mynu] = SpmCorrMatZERO(eigenvalsp, eigenvecsp, converttablep, beta, time);
      stop = clock();

      if (Nh > 1)
      {
        corrfuncNhZERO[mynu] = NhCorrZERO(beta, time);
      }

      cout << "SpmMAT time for nu = " << mynu << ": " << (stop-start)/CLOCKS_PER_SEC << endl;
    }
    if(NNCORR)
    {
      Sminus = makeSminus(converttablep);

      cout << "BEFORE SZCORRNN FOR nu=" << nu << endl;

      start = clock(); //start clock
      corrfunczzNN[mynu] = SzCorrMatNN(beta, time);
      stop = clock();

      cout << "SzMATNN time for nu = " << mynu << ": " << (stop-start)/CLOCKS_PER_SEC << endl;

      start = clock(); //start clock
      corrfuncpmNN[mynu] = SpmCorrMatNN(eigenvalsp, eigenvecsp, converttablep, beta, time);
      stop = clock();

      if (Nh > 1)
      {
        corrfuncNhNN[mynu] = NhCorrNN(beta, time);
      }

      cout << "SpmMATNN time for nu = " << mynu << ": " << (stop-start)/CLOCKS_PER_SEC << endl;
    }
    if(MIDCORR)
    {
      Sminus = makeSminus(converttablep);

      cout << "BEFORE SZCORRMID FOR nu=" << nu << endl;

      start = clock(); //start clock
      corrfunczzMID[mynu] = SzCorrMatMID(beta, time);
      stop = clock();

      cout << "SzMATMID time for nu = " << mynu << ": " << (stop-start)/CLOCKS_PER_SEC << endl;

      start = clock(); //start clock
      corrfuncpmMID[mynu] = SpmCorrMatMID(eigenvalsp, eigenvecsp, converttablep, beta, time);
      stop = clock();

      if (Nh > 1)
      {
        corrfuncNhMID[mynu] = NhCorrMID(beta, time);
      }

      cout << "SpmMATMID time for nu = " << mynu << ": " << (stop-start)/CLOCKS_PER_SEC << endl;
    }
    eigenvalsp = eigenvals;
    eigenvecsp = eigenvecs;
    converttablep = converttable;
  }

  double GS = findminimum(mineigvals);
  cout << GS << endl;
  vector<double> partitionfunction(Nb, 0.0);

  vector<vector<vector<complex<double>>>> corrzZERO(TWOL, vector<vector<complex<double>>>(Nb, vector<complex<double>>(Nt, zero)));
  vector<vector<vector<complex<double>>>> corrpmZERO(TWOL, vector<vector<complex<double>>>(Nb, vector<complex<double>>(Nt, zero)));
  vector<vector<vector<complex<double>>>> corrNhZERO(TWOL, vector<vector<complex<double>>>(Nb, vector<complex<double>>(Nt, zero)));

  vector<vector<vector<complex<double>>>> corrzNN(TWOL, vector<vector<complex<double>>>(Nb, vector<complex<double>>(Nt, zero)));
  vector<vector<vector<complex<double>>>> corrpmNN(TWOL, vector<vector<complex<double>>>(Nb, vector<complex<double>>(Nt, zero)));
  vector<vector<vector<complex<double>>>> corrNhNN(TWOL, vector<vector<complex<double>>>(Nb, vector<complex<double>>(Nt, zero)));

  vector<vector<vector<complex<double>>>> corrzMID(TWOL, vector<vector<complex<double>>>(Nb, vector<complex<double>>(Nt, zero)));
  vector<vector<vector<complex<double>>>> corrpmMID(TWOL, vector<vector<complex<double>>>(Nb, vector<complex<double>>(Nt, zero)));
  vector<vector<vector<complex<double>>>> corrNhMID(TWOL, vector<vector<complex<double>>>(Nb, vector<complex<double>>(Nt, zero)));

  cout << endl;
  cout << endl;

  for(int i = 0; i < Ns+1; i++)
    for(int b = 0; b < Nb; b++)
    {
      //cout << partfunc[i]*exp(-beta*(mineigvals[i]-GS)) << endl;
      partitionfunction[b] += partfunc[i][b]*exp(-beta[b]*(mineigvals[i]-GS));//Sould double check this for a small system?
      for(int t = 0; t < Nt; t++)
        for(int j = 0; j < TWOL; j++)
        {
          //cout << "Here: " << corrfuncpm[i][j] << "   " << exp(-beta*(mineigvals[i]-GS)) << endl;
          corrzZERO[j][b][t] += corrfunczzZERO[i][j][b][t]*exp(-beta[b]*(mineigvals[i]-GS));
          corrpmZERO[j][b][t] += corrfuncpmZERO[i][j][b][t]*exp(-beta[b]*(mineigvalspm[i]-GS));
          corrNhZERO[j][b][t] += corrfuncNhZERO[i][j][b][t]*exp(-beta[b]*(mineigvals[i]-GS));
          corrzNN[j][b][t] += corrfunczzNN[i][j][b][t]*exp(-beta[b]*(mineigvals[i]-GS));
          corrpmNN[j][b][t] += corrfuncpmNN[i][j][b][t]*exp(-beta[b]*(mineigvalspm[i]-GS));
          corrNhNN[j][b][t] += corrfuncNhNN[i][j][b][t]*exp(-beta[b]*(mineigvals[i]-GS));
          corrzMID[j][b][t] += corrfunczzMID[i][j][b][t]*exp(-beta[b]*(mineigvals[i]-GS));
          corrpmMID[j][b][t] += corrfuncpmMID[i][j][b][t]*exp(-beta[b]*(mineigvalspm[i]-GS));
          corrNhMID[j][b][t] += corrfuncNhMID[i][j][b][t]*exp(-beta[b]*(mineigvals[i]-GS));
        }
    }

  for(int j = 0; j < TWOL; j++)
    for(int b = 0; b < Nb; b++)
      for(int t = 0; t < Nt; t++)
      {
        //cout << "Here: " << corrfuncpm[i][j] << "   " << exp(-beta*(mineigvals[i]-GS)) << endl;
        corrzZERO[j][b][t] /= partitionfunction[b];
        corrpmZERO[j][b][t] /= partitionfunction[b];
        corrNhZERO[j][b][t] /= partitionfunction[b];
        corrzNN[j][b][t] /= partitionfunction[b];
        corrpmNN[j][b][t] /= partitionfunction[b];
        corrNhNN[j][b][t] /= partitionfunction[b];
        corrzMID[j][b][t] /= partitionfunction[b];
        corrpmMID[j][b][t] /= partitionfunction[b];
        corrNhMID[j][b][t] /= partitionfunction[b];
      }

  cout << "Feridg" << endl;

  if(ZEROCORR)
  {
    WriteCorr(beta, time, corrzZERO, corrpmZERO, "Corr.txt");
    if (Nh > 1) WriteHoleCorr(beta, time, corrNhZERO, "HoleCorr.txt");
  }
  if(NNCORR)
  {
    WriteCorr(beta, time, corrzNN, corrpmNN, "CorrNN.txt");
    if (Nh > 1) WriteHoleCorr(beta, time, corrNhNN, "HoleCorrNN.txt");
  }
  if(MIDCORR)
  {
    WriteCorr(beta, time, corrzMID, corrpmMID, "CorrMID.txt");
    if (Nh > 1) WriteHoleCorr(beta, time, corrNhMID, "HoleCorrMID.txt");
  }

  WritePartition(beta, partitionfunction);

  cout << "Ferdig" << endl;
}


void Solver::makebasis() //nu is number of up spins
{
  maxIndexValue = binomial(TWOL, Nh)*binomial(Ns, nu); //This is the number of states with Nh holes and nu up-spins.
  converttable.init_index_to_state(maxIndexValue);

  //Construct all states with nu up spins.

  //How do I want to represent a state? Possibly as a vector with 2L numbers -1, 0 or 1? Or as a vector of length Nh and a binary number with the rest?
  //Do it as a vector of length Nh and a binary number with the rest and translate it into spins on sites when computing quantities.
  //For 0 and 1 hole I can do s + 2^(2L-Nh)*holeposition. But what to do for two or more holes? Think about this for later.

  //First find all the possible binary numbers with nu ones. Loop over all numbers from 0 to 2^{2*L-Nh}

  vector<unsigned int> s = {}; //List of all spin configurations with nu up spins.

  //int nconf = 0;
  for (int i = 0; i < twomax; i++)
  {
    int nunow = 0;
    for (int j = 0; j < Ns; j++)
    {
      if (i >> j & 1) nunow += 1;
    }
    if (nunow == nu)
    {
      //cout << dectobin(i) << endl;
      //nconf += 1;
      s.push_back(i);
    }
  }

  int NsStates = s.size(); //Seems to be correct.

  short unsigned int index;
  unsigned int state;
  int maxStateValue = 0;

  //Then combine these with all possible hole positions to create a number representing a state. And make lists converting between state numbers and matrix indices.
  //Make an algorithm for general Nh:

  vector<int> holeind(Nh+1);
  vector<int> holemax(Nh+1);

  for (int i = 0; i < Nh; i++)
  {
    holeind[i] = Nh-1-i;
    holemax[i] = TWOL-i;
  }
  holeind[Nh] = 0;
  holemax[Nh] = 2;

  int p = 0;

  index = 0;
  while (holeind[Nh] == 0)
  {
    /*for (int kk = 0; kk < Nh; kk++)
    {
      cout << holeind[kk] << "   ";
    }
    cout << endl;*/
    for (int i = 0; i < NsStates; i++)
    {
      state = s[i];
      for (int j = 0; j < Nh; j++)
      {
        state += twomax*(TWOLpow[j]*holeind[Nh-1-j]); //NOT SURE WHETHER I AM COMBINING THE CORRECT TWOLpow AND ind?!
        //cout << TWOLpow[j] << "   " << holeind[Nh-1-j] << endl;
      }
      if (state > maxStateValue) maxStateValue = state;
      converttable.index_to_state[index] = state;
      //cout << state << endl;
      index++;
    }

    holeind[0]++;
    while (holeind[p] == holemax[p])
    {
      ++p;
      //cout << "p = " << p << endl;
      holeind[p]++;
      for (int q = p-1; q>=0; q--) holeind[q]=holeind[q+1]+1; //Feilen ligger her


      if(holeind[Nh] == 1){break;}
      if(holeind[p]!=holemax[p]){p=0;}
    }
  }

  converttable.init_state_to_index(maxStateValue+1);
  for (int i = 0; i < maxStateValue; i++) converttable.state_to_index[i] = maxIndexValue;
  for (int i = 0; i < maxIndexValue; i++)
  {
    converttable.state_to_index[converttable.index_to_state[i]] = i;
  }

  /*
  for (int i = 0; i < converttable.state_to_index.size(); i++)
  {
    cout << converttable.state_to_index[i] << endl;
  }

  for (int i = 0; i < converttable.index_to_state.size(); i++)
  {
    cout << converttable.index_to_state[i] << endl;
  }
  */

  //I think these work??? Look over it again tomorrow.


  //What to do about all the zeros which are not supposed to be zeros? I have set them to the maxIndexValue, i.e. out of range.
}

unsigned int Solver::statevec_to_statenum(vector<short int> statevec)
{
  int statenum = 0;

  int spin = 0;
  int holeexp = 0;

  for (int i = 0; i < statevec.size(); i++)
  {
    if (statevec[i] == 0)
    {
      statenum += i*twomax*pow(TWOL, holeexp);
      holeexp++;
      continue;
    }
    spin = spin << 1;
    if (statevec[i] == 1)
    {
      spin += 1;
    }
  }
  statenum += spin;
  return statenum;
}

vector<short int> Solver::statenum_to_statevec(unsigned int statenum)
{
  vector<short int> statevec(TWOL, -2);

  unsigned int spin = statenum%twomax;

  statenum /= twomax;
  //cout << statenum << endl;

  for (int i = 0; i < Nh; i++)
  {
    statevec[statenum/TWOLpow[Nh-1-i]] = 0;
    statenum %= TWOLpow[Nh-1-i];
  }

  for (int i = TWOL-1; i >= 0; i--)
  {
    if (statevec[i] != 0)
    {
      statevec[i] += 1 + 2*(spin%2);
      spin = spin >> 1;
    }
  }

  return statevec;
}

void Solver::fillH()
{
  //Compute all elements of H for a given nu?

  //maxIndexValue = binomial(2*L, Nh)*binomial(Ns, nu);

  H.resize(maxIndexValue,maxIndexValue);
  H.setZero();

  unsigned int innumber;
  vector<short int> instate;
  vector<short int> holevec(Nh);
  int holepos;
  int neighpos;
  int neigh2;
  int additionalsign;

  int j;

  for (int i = 0; i < maxIndexValue; i++)
  {

    innumber = converttable.index_to_state[i];
    instate = statenum_to_statevec(innumber);

    innumber /= twomax;

    for (int h = 0; h < Nh; h++)
    {
      holevec[Nh-1-h] = innumber/TWOLpow[Nh-1-h];
      innumber %= TWOLpow[Nh-1-h];
    }

    if (PYROCHLORE)
    {
      //Diagonal terms are easy: SziSzj for Jleg and Jrung add a Delta to tune z-component?:

      for (int site = 0; site < TWOL-1; site++)
      {
        H(i, i) += 0.25*Jzr*instate[site]*instate[site+1];
        if (site%3 == 0) H(i, i) += 0.25*Jzr*instate[site]*instate[site+3];
      }
      for (int site = 0; site < TWOL-2; site++)
      {
        if ((site-2)%3 != 0) H(i, i) += 0.25*Jzl*instate[site]*instate[site+2];
      }


      //Off-diagonal terms:

      //Heisenberg part:
      //Nearest neighbours may flip if they are one up and one down.

      for (int site = 0; site < TWOL-1; site++)
      {
        if (instate[site]*instate[site+1] == -1)
        {
          instate[site] *= (-1);
          instate[site+1] *= (-1);

          j = converttable.state_to_index[statevec_to_statenum(instate)];

          H(i,j) += 0.5*Jpmr;

          instate[site] *= (-1);
          instate[site+1] *= (-1);
        }
      }
      for (int site = 0; site < TWOL-3; site+=3)
      {
        if (instate[site]*instate[site+3] == -1)
        {
          instate[site] *= (-1);
          instate[site+3] *= (-1);

          j = converttable.state_to_index[statevec_to_statenum(instate)];

          H(i,j) += 0.5*Jpmr;

          instate[site] *= (-1);
          instate[site+3] *= (-1);
        }
      }
      for (int site = 0; site < TWOL-2; site++)
      {
        if (instate[site]*instate[site+2] == -1)
        {
          instate[site] *= (-1);
          instate[site+2] *= (-1);

          j = converttable.state_to_index[statevec_to_statenum(instate)];

          if ((site-2)%3 != 0) H(i,j) += 0.5*Jpml;

          instate[site] *= (-1);
          instate[site+2] *= (-1);
        }
      }

      //Hopping part:
      for (int neigh = -1; neigh < 2; neigh += 2)
      {
        neigh2 = 2*neigh;
        for (int hole = 0; hole < Nh; hole++)
        {
          holepos = holevec[hole];

          //Hop along rungs
          neighpos = (holepos+neigh+TWOL)%TWOL;

          if ((neigh == -1 && holepos == 0) || (neigh == +1 && holepos == TWOL-1)) {additionalsign = 0;}
          else {additionalsign = +1;}

          instate[holepos] = instate[neighpos];
          instate[neighpos] = 0;
          j = converttable.state_to_index[statevec_to_statenum(instate)];

          H(i,j) += -tr*additionalsign*abs(instate[holepos]);

          instate[neighpos] = instate[holepos];
          instate[holepos] = 0;


          //Hop along legs.
          neighpos = (holepos+neigh2+TWOL)%TWOL;

          if ((neigh2 == -2 && holepos == 0) || (neigh2 == +2 && holepos == TWOL-2) || (neigh2 == -2 && holepos == 1) || (neigh2 == +2 && holepos == TWOL-1)) {additionalsign = 0;}
          else if ((neigh2 == +2 && (holepos-2)%3 == 0) || (neigh2 == -2 && (holepos-4)%3 == 0)) {additionalsign = 0;}
          else if (instate[(holepos+neigh+TWOL)%TWOL] == 0) {additionalsign = +1;}
          else {additionalsign = -1;}

          instate[holepos] = instate[neighpos];
          instate[neighpos] = 0;
          j = converttable.state_to_index[statevec_to_statenum(instate)];
          H(i,j) += -tl*additionalsign*abs(instate[holepos]);
          //cout << i << "   " << j << "   " << -tl*additionalsign*abs(instate[holepos]) << endl;

          instate[neighpos] = instate[holepos];
          instate[holepos] = 0;

          //Hop along third edge of tetrahedron.
          if (holepos%3 == 0)
          {
            int neigh3 = 3*neigh;
            neighpos = (holepos+neigh3+TWOL)%TWOL;

            //cout << holepos << "    " << neighpos << endl;

            if ((neigh3 == -3 && holepos == 0) || (neigh3 == +3 && holepos == TWOL-1)) {additionalsign = 0;}
            else {additionalsign = pow((-1),(abs(instate[holepos+1*neigh])+abs(instate[holepos+2*neigh])));}

            //cout << additionalsign << endl;

            instate[holepos] = instate[neighpos];
            instate[neighpos] = 0;
            j = converttable.state_to_index[statevec_to_statenum(instate)];
            H(i,j) += -tr*additionalsign*abs(instate[holepos]);
            //cout << i << "   " << j << "   " << -tl*additionalsign*abs(instate[holepos]) << endl;

            instate[neighpos] = instate[holepos];
            instate[holepos] = 0;
          }
        }
      }
    }
    else if(OBC)
    {
      //Diagonal terms are easy: SziSzj for Jleg and Jrung add a Delta to tune z-component?:

      for (int site = 0; site < TWOL-1; site++)
      {
        H(i, i) += 0.25*Jzr*instate[site]*instate[site+1];
      }
      for (int site = 0; site < TWOL-2; site++)
      {
        H(i, i) += 0.25*Jzl*instate[site]*instate[site+2];
      }


      //Off-diagonal terms:

      //Heisenberg part:
      //Nearest neighbours may flip if they are one up and one down.

      for (int site = 0; site < TWOL-1; site++)
      {
        if (instate[site]*instate[site+1] == -1)
        {
          instate[site] *= (-1);
          instate[site+1] *= (-1);

          j = converttable.state_to_index[statevec_to_statenum(instate)];

          H(i,j) += 0.5*Jpmr;

          instate[site] *= (-1);
          instate[site+1] *= (-1);
        }
      }
      for (int site = 0; site < TWOL-2; site++)
      {
        if (instate[site]*instate[site+2] == -1)
        {
          instate[site] *= (-1);
          instate[site+2] *= (-1);

          j = converttable.state_to_index[statevec_to_statenum(instate)];

          H(i,j) += 0.5*Jpml;

          instate[site] *= (-1);
          instate[site+2] *= (-1);
        }
      }

      //Hopping part:
      for (int neigh = -1; neigh < 2; neigh += 2)
      {
        neigh2 = 2*neigh;
        for (int hole = 0; hole < Nh; hole++)
        {
          holepos = holevec[hole];

          //Hop along rungs
          neighpos = (holepos+neigh+TWOL)%TWOL;

          if ((neigh == -1 && holepos == 0) || (neigh == +1 && holepos == TWOL-1)) {additionalsign = 0;}
          else {additionalsign = +1;}

          instate[holepos] = instate[neighpos];
          instate[neighpos] = 0;
          j = converttable.state_to_index[statevec_to_statenum(instate)];

          H(i,j) += -tr*additionalsign*abs(instate[holepos]);

          instate[neighpos] = instate[holepos];
          instate[holepos] = 0;


          //Hop along legs.
          neighpos = (holepos+neigh2+TWOL)%TWOL;

          if ((neigh2 == -2 && holepos == 0) || (neigh2 == +2 && holepos == TWOL-2) || (neigh2 == -2 && holepos == 1) || (neigh2 == +2 && holepos == TWOL-1)) {additionalsign = 0;}
          else if (instate[(holepos+neigh+TWOL)%TWOL] == 0) {additionalsign = +1;}
          else {additionalsign = -1;}

          instate[holepos] = instate[neighpos];
          instate[neighpos] = 0;
          j = converttable.state_to_index[statevec_to_statenum(instate)];
          H(i,j) += -tl*additionalsign*abs(instate[holepos]);
          //cout << i << "   " << j << "   " << -tl*additionalsign*abs(instate[holepos]) << endl;

          instate[neighpos] = instate[holepos];
          instate[holepos] = 0;
        }
      }
    }
    else
    {
      //Diagonal terms are easy: SziSzj for Jleg and Jrung add a Delta to tune z-component?:

      for (int site = 0; site < TWOL; site++)
      {
        H(i, i) += 0.25*2*Jzr*instate[site]*instate[(site+1)%TWOL];
        H(i, i) += 0.25*Jzl*instate[site]*instate[(site+2)%TWOL]; //Should I include the 0.25 from s=1/2?
      }


      //Off-diagonal terms:

      //Heisenberg part:
      //Nearest neighbours may flip if they are one up and one down.

      for (int site = 0; site < TWOL; site++)
      {
        if (instate[site]*instate[(site+1)%TWOL] == -1)
        {
          instate[site] *= (-1);
          instate[(site+1)%TWOL] *= (-1);

          j = converttable.state_to_index[statevec_to_statenum(instate)];

          H(i,j) += 0.5*2*Jpmr;

          instate[site] *= (-1);
          instate[(site+1)%TWOL] *= (-1);
        }
        if (instate[site]*instate[(site+2)%TWOL] == -1)
        {
          instate[site] *= (-1);
          instate[(site+2)%TWOL] *= (-1);

          j = converttable.state_to_index[statevec_to_statenum(instate)];

          H(i,j) += 0.5*Jpml;

          instate[site] *= (-1);
          instate[(site+2)%TWOL] *= (-1);
        }
      }

      //Hopping part:
      for (int neigh = -1; neigh < 2; neigh += 2)
      {
        neigh2 = 2*neigh;
        for (int hole = 0; hole < Nh; hole++)
        {
          holepos = holevec[hole];

          //Hop along rungs
          neighpos = (holepos+neigh+TWOL)%TWOL;

          if ((neigh == -1 && holepos == 0) || (neigh == +1 && holepos == TWOL-1)) {additionalsign = +1 + 2*(-1)*((TWOL-2-(Nh-1))%2);}
          else {additionalsign = +1;}

          instate[holepos] = instate[neighpos];
          instate[neighpos] = 0;
          j = converttable.state_to_index[statevec_to_statenum(instate)];

          H(i,j) += -2*tr*additionalsign*abs(instate[holepos]);

          instate[neighpos] = instate[holepos];
          instate[holepos] = 0;


          //Hop along legs.
          neighpos = (holepos+neigh2+TWOL)%TWOL;

          if ((neigh2 == -2 && holepos == 0) || (neigh2 == +2 && holepos == TWOL-2)){additionalsign = (+1) + 2*(-1)*((TWOL-3-(Nh-1+(int(abs(instate[TWOL-1]))-1)))%2);}
          else if ((neigh2 == -2 && holepos == 1) || (neigh2 == +2 && holepos == TWOL-1)) {additionalsign = (+1) + 2*(-1)*((TWOL-3-(Nh-1+(int(abs(instate[0]))-1)))%2);}
          else if (instate[(holepos+neigh+TWOL)%TWOL] == 0) {additionalsign = +1;}
          else {additionalsign = -1;}

          instate[holepos] = instate[neighpos];
          instate[neighpos] = 0;
          j = converttable.state_to_index[statevec_to_statenum(instate)];
          H(i,j) += -tl*additionalsign*abs(instate[holepos]);
          //cout << i << "   " << j << "   " << -tl*additionalsign*abs(instate[holepos]) << endl;

          instate[neighpos] = instate[holepos];
          instate[holepos] = 0;
        }
      }
    }
  }


  return;
}

void Solver::diagonalise()
{
  SelfAdjointEigenSolver<Matrix<double,Dynamic,Dynamic>> es;
  if (EIGVEC)
  {
    es = SelfAdjointEigenSolver<Matrix<double,Dynamic,Dynamic>>(H,ComputeEigenvectors);
    eigenvecs = es.eigenvectors();
  }
  else {es = SelfAdjointEigenSolver<Matrix<double,Dynamic,Dynamic>>(H,EigenvaluesOnly);}

  eigenvals = es.eigenvalues();

  return;
}

Matrix<double,Dynamic,Dynamic> Solver::NhMATRIX(int siteind)
{
  Matrix<double,Dynamic,Dynamic> ans(maxIndexValue,maxIndexValue);
  double statenumber;
  vector<short int> statevec;

  ans.setZero();

  for (int i = 0; i < maxIndexValue; i++)
  {
    statenumber = converttable.index_to_state[i];
    statevec = statenum_to_statevec(statenumber);
    ans(i,i) = 1-abs(statevec[siteind]);
  }

  return ans;
}

vector<vector<complex<double>>> Solver::NhExp(vector<double> beta, vector<double> time)
{
  vector<Matrix<double,Dynamic,Dynamic>> Nhjeigenbasis;

  Matrix<double,Dynamic,Dynamic> eigenvecsinv;

  eigenvecsinv = eigenvecs.inverse();

  for(int j = 0; j < TWOL; j++)
  {
    Nhjeigenbasis.push_back(eigenvecsinv*NhMATRIX(j)*eigenvecs);
  }

  double minval = eigenvals[0];

  vector<vector<complex<double>>> msum(TWOL, vector<complex<double>>(Nt, zero));
  vector<vector<complex<double>>> NhExpVal(TWOL, vector<complex<double>>(Nb, zero));

  for(int n = 0; n < maxIndexValue; n++)
    for(int j = 0; j < TWOL; j++)
      for(int b = 0; b < Nb; b++)
      {
        NhExpVal[j][b] += exp(-beta[b]*(eigenvals[n]-minval))*Nhjeigenbasis[j](n, n);
      }

  return NhExpVal;
}

vector<vector<vector<complex<double>>>> Solver::NhCorrZERO(vector<double> beta, vector<double> time)
{
  Matrix<double,Dynamic,Dynamic> Nh0eigenbasis;
  vector<Matrix<double,Dynamic,Dynamic>> Nhjeigenbasis;

  Matrix<double,Dynamic,Dynamic> eigenvecsinv;

  eigenvecsinv = eigenvecs.inverse();

  for(int j = 0; j < TWOL; j++)
  {
    Nhjeigenbasis.push_back(eigenvecsinv*NhMATRIX(j)*eigenvecs);
  }
  Nh0eigenbasis = Nhjeigenbasis[0];

  double minval = eigenvals[0];

  vector<vector<complex<double>>> msum(TWOL, vector<complex<double>>(Nt, zero));
  vector<vector<vector<complex<double>>>> ans(TWOL, vector<vector<complex<double>>>(Nb, vector<complex<double>>(Nt, zero)));

  for(int n = 0; n < maxIndexValue; n++)
  {
    for(int j = 0; j < TWOL; j++)
      for(int t = 0; t < Nt; t++)
      {
        msum[j][t] = zero;
        for(int m = 0; m < maxIndexValue; m++)
        {
          msum[j][t] += Nh0eigenbasis(n, m)*Nhjeigenbasis[j](m, n)*exponential(-(eigenvals[n]-eigenvals[m])*time[t]);
        }
      }
    for(int j = 0; j < TWOL; j++)
      for(int t = 0; t < Nt; t++)
        for(int b = 0; b < Nb; b++)
        {
          ans[j][b][t] += exp(-beta[b]*(eigenvals[n]-minval))*msum[j][t];
        }
  }
  return ans;
}

vector<vector<vector<complex<double>>>> Solver::NhCorrNN(vector<double> beta, vector<double> time)
{
  vector<Matrix<double,Dynamic,Dynamic>> Nhjeigenbasis;

  Matrix<double,Dynamic,Dynamic> eigenvecsinv;

  eigenvecsinv = eigenvecs.inverse();

  for(int j = 0; j < TWOL; j++)
  {
    Nhjeigenbasis.push_back(eigenvecsinv*NhMATRIX(j)*eigenvecs);
  }

  double minval = eigenvals[0];

  vector<vector<complex<double>>> msum(TWOL, vector<complex<double>>(Nt, zero));
  vector<vector<vector<complex<double>>>> ans(TWOL, vector<vector<complex<double>>>(Nb, vector<complex<double>>(Nt, zero)));

  for(int n = 0; n < maxIndexValue; n++)
  {
    for(int j = 0; j < TWOL; j++)
      for(int t = 0; t < Nt; t++)
      {
        msum[j][t] = zero;
        for(int m = 0; m < maxIndexValue; m++)
        {
          msum[j][t] += Nhjeigenbasis[j](n, m)*Nhjeigenbasis[(j+1)%TWOL](m, n)*exponential(-(eigenvals[n]-eigenvals[m])*time[t]);
        }
      }
    for(int j = 0; j < TWOL; j++)
      for(int t = 0; t < Nt; t++)
        for(int b = 0; b < Nb; b++)
        {
          ans[j][b][t] += exp(-beta[b]*(eigenvals[n]-minval))*msum[j][t];
        }
  }
  return ans;
}

vector<vector<vector<complex<double>>>> Solver::NhCorrMID(vector<double> beta, vector<double> time)
{
  Matrix<double,Dynamic,Dynamic> NhMIDeigenbasis;
  vector<Matrix<double,Dynamic,Dynamic>> Nhjeigenbasis;

  Matrix<double,Dynamic,Dynamic> eigenvecsinv;

  eigenvecsinv = eigenvecs.inverse();

  for(int j = 0; j < TWOL; j++)
  {
    Nhjeigenbasis.push_back(eigenvecsinv*NhMATRIX(j)*eigenvecs);
  }
  NhMIDeigenbasis = Nhjeigenbasis[int(Nsites/2)];

  double minval = eigenvals[0];

  vector<vector<complex<double>>> msum(TWOL, vector<complex<double>>(Nt, zero));
  vector<vector<vector<complex<double>>>> ans(TWOL, vector<vector<complex<double>>>(Nb, vector<complex<double>>(Nt, zero)));

  for(int n = 0; n < maxIndexValue; n++)
  {
    for(int j = 0; j < TWOL; j++)
      for(int t = 0; t < Nt; t++)
      {
        msum[j][t] = zero;
        for(int m = 0; m < maxIndexValue; m++)
        {
          msum[j][t] += NhMIDeigenbasis(n, m)*Nhjeigenbasis[j](m, n)*exponential(-(eigenvals[n]-eigenvals[m])*time[t]);
        }
      }
    for(int j = 0; j < TWOL; j++)
      for(int t = 0; t < Nt; t++)
        for(int b = 0; b < Nb; b++)
        {
          ans[j][b][t] += exp(-beta[b]*(eigenvals[n]-minval))*msum[j][t];
        }
  }
  return ans;
}

vector<double> Solver::Szmat(int siteind)
{
  //Sovle it as a sparse matrix maybe? It is diagonal?
  vector<double> ans(maxIndexValue);
  double statenumber;
  vector<short int> statevec;

  for (int i = 0; i < maxIndexValue; i++)
  {
    statenumber = converttable.index_to_state[i];
    statevec = statenum_to_statevec(statenumber);
    ans[i] = 0.5*statevec[siteind];
  }

  return ans;
}

Matrix<double,Dynamic,Dynamic> Solver::SzmatMATRIX(int siteind)
{
  //Sovle it as a sparse matrix maybe? It is diagonal?
  Matrix<double,Dynamic,Dynamic> ans(maxIndexValue,maxIndexValue);
  double statenumber;
  vector<short int> statevec;

  ans.setZero();

  for (int i = 0; i < maxIndexValue; i++)
  {
    statenumber = converttable.index_to_state[i];
    statevec = statenum_to_statevec(statenumber);
    ans(i,i) = 0.5*statevec[siteind];
  }

  return ans;
}


vector<vector<vector<complex<double>>>> Solver::SzCorrMatZERO(vector<double> beta, vector<double> time)
{
  Matrix<double,Dynamic,Dynamic> S0eigenbasis;
  vector<Matrix<double,Dynamic,Dynamic>> Sjeigenbasis;

  Matrix<double,Dynamic,Dynamic> eigenvecsinv;

  double start2 = clock();

  eigenvecsinv = eigenvecs.inverse();

  double stop2 = clock();

  cout << "SzMAT invert time: " << (stop2-start2)/CLOCKS_PER_SEC << endl;

  for(int j = 0; j < TWOL; j++)
  {
    Sjeigenbasis.push_back(eigenvecsinv*SzmatMATRIX(j)*eigenvecs);
  }
  S0eigenbasis = Sjeigenbasis[0];

  double minval = eigenvals[0];

  vector<vector<complex<double>>> msum(TWOL, vector<complex<double>>(Nt, zero));
  vector<vector<vector<complex<double>>>> SzSz(TWOL, vector<vector<complex<double>>>(Nb, vector<complex<double>>(Nt, zero)));

  double start1 = clock(); //start clock
  if(CUTOFF)
  {
    //Compute correlations only for energies up to T bigger than GS.
    //Doing only static to save time.
    int n;

    for(int b = 0; b < Nb; b++)
    {
      n = 0;
      while(abs(eigenvals[n]-minval) < 1./beta[b]) //Only for n-sum?
      {
        for(int j = 0; j < TWOL; j++)
          for(int t = 0; t < Nt; t++)
          {
            msum[j][t] = zero;
            for(int m = 0; m < maxIndexValue; m++)
            {
              msum[j][t] += S0eigenbasis(n, m)*Sjeigenbasis[j](m, n);
            }
          }
        for(int j = 0; j < TWOL; j++)
          for(int t = 0; t < Nt; t++)
          {
            SzSz[j][b][t] += exp(-beta[b]*(eigenvals[n]-minval))*msum[j][t];
          }
        n++;
      }
    }

  }
  else
  {
    for(int n = 0; n < maxIndexValue; n++) //Only for n-sum?
    {
      for(int j = 0; j < TWOL; j++)
        for(int t = 0; t < Nt; t++)
        {
          msum[j][t] = zero;
          for(int m = 0; m < maxIndexValue; m++)
          {
            msum[j][t] += S0eigenbasis(n, m)*Sjeigenbasis[j](m, n)*exponential(-(eigenvals[n]-eigenvals[m])*time[t]);
          }
        }
      for(int j = 0; j < TWOL; j++)
        for(int t = 0; t < Nt; t++)
          for(int b = 0; b < Nb; b++)
          {
            SzSz[j][b][t] += exp(-beta[b]*(eigenvals[n]-minval))*msum[j][t];
          }
    }
  }

  double stop1 = clock();

  cout << "SzMAT core time: " << (stop1-start1)/CLOCKS_PER_SEC << endl;

  return SzSz;
}

vector<vector<vector<complex<double>>>> Solver::SzCorrMatNN(vector<double> beta, vector<double> time)
{
  vector<Matrix<double,Dynamic,Dynamic>> Sjeigenbasis;

  Matrix<double,Dynamic,Dynamic> eigenvecsinv;

  eigenvecsinv = eigenvecs.inverse();

  for(int j = 0; j < TWOL; j++)
  {
    Sjeigenbasis.push_back(eigenvecsinv*SzmatMATRIX(j)*eigenvecs);
  }

  double minval = eigenvals[0];

  vector<vector<complex<double>>> msum(TWOL, vector<complex<double>>(Nt, zero));
  vector<vector<vector<complex<double>>>> SzSz(TWOL, vector<vector<complex<double>>>(Nb, vector<complex<double>>(Nt, zero)));

  int maxInd;
  if(CUTOFF){maxInd = maxIndexValue;}
  else {maxInd = maxIndexValue;}

  for(int n = 0; n < maxInd; n++) //Only for n-sum?
  {
    for(int j = 0; j < TWOL; j++)
      for(int t = 0; t < Nt; t++)
      {
        msum[j][t] = zero;
        for(int m = 0; m < maxIndexValue; m++)
        {
          msum[j][t] += Sjeigenbasis[j](n, m)*Sjeigenbasis[(j+1)%TWOL](m, n)*exponential(-(eigenvals[n]-eigenvals[m])*time[t]);
        }
      }
    for(int j = 0; j < TWOL; j++)
      for(int t = 0; t < Nt; t++)
        for(int b = 0; b < Nb; b++)
        {
          SzSz[j][b][t] += exp(-beta[b]*(eigenvals[n]-minval))*msum[j][t];
        }
  }
  return SzSz;
}

vector<vector<vector<complex<double>>>> Solver::SzCorrMatMID(vector<double> beta, vector<double> time)
{
  Matrix<double,Dynamic,Dynamic> SMIDeigenbasis;
  vector<Matrix<double,Dynamic,Dynamic>> Sjeigenbasis;

  Matrix<double,Dynamic,Dynamic> eigenvecsinv;

  eigenvecsinv = eigenvecs.inverse();

  for(int j = 0; j < TWOL; j++)
  {
    Sjeigenbasis.push_back(eigenvecsinv*SzmatMATRIX(j)*eigenvecs);
  }
  SMIDeigenbasis = Sjeigenbasis[int(Nsites/2)];

  double minval = eigenvals[0];

  vector<vector<complex<double>>> msum(TWOL, vector<complex<double>>(Nt, zero));
  vector<vector<vector<complex<double>>>> SzSz(TWOL, vector<vector<complex<double>>>(Nb, vector<complex<double>>(Nt, zero)));

  int maxInd;
  if(CUTOFF){maxInd = maxIndexValue;}
  else {maxInd = maxIndexValue;}

  for(int n = 0; n < maxInd; n++) //Only for n-sum?
  {
    for(int j = 0; j < TWOL; j++)
      for(int t = 0; t < Nt; t++)
      {
        msum[j][t] = zero;
        for(int m = 0; m < maxIndexValue; m++)
        {
          msum[j][t] += SMIDeigenbasis(n, m)*Sjeigenbasis[j](m, n)*exponential(-(eigenvals[n]-eigenvals[m])*time[t]);
        }
      }
    for(int j = 0; j < TWOL; j++)
      for(int t = 0; t < Nt; t++)
        for(int b = 0; b < Nb; b++)
        {
          SzSz[j][b][t] += exp(-beta[b]*(eigenvals[n]-minval))*msum[j][t];
        }
  }
  return SzSz;
}

vector<vector<vector<complex<double>>>> Solver::SzCorr(vector<double> beta, vector<double> time)
{
  //vector<vector<complex<double>>> W(Nb, vector<complex<double>>(Nt, zero));
  //vector<vector<vector<double>>>partsum(TWOL, vector<vector<double>>(maxIndexValue, vector<double>(maxIndexValue, 0.0)));
  vector<vector<vector<complex<double>>>> SzSz(TWOL, vector<vector<complex<double>>>(Nb, vector<complex<double>>(Nt, zero)));

  double nSz0m = 0;
  double mSzjn = 0;

  double minval = eigenvals[0];

  //Compute SzSz correlations between site i and all other sites.
  vector<double> Sz0 = Szmat(0);
  vector<vector<double>> Szj(TWOL, vector<double>(maxIndexValue, 0.0));
  for(int j = 0; j < TWOL; j++) Szj[j] = Szmat(j);

  Eigen::Matrix<double, -1, 1, 0, -1, 1> eigenvecn, eigenvecm;
  //double eigenvecnmalpha;
  vector<double> Szjnow(maxIndexValue);

  for (int n = 0; n < maxIndexValue; n++)
  {
    eigenvecn = eigenvecs.col(n);
    for (int m = 0; m < maxIndexValue; m++)
    {
      eigenvecm = eigenvecs.col(m);
      nSz0m = 0;
      for (int alpha = 0; alpha < maxIndexValue; alpha ++)
      {
        nSz0m += eigenvecn[alpha]*eigenvecm[alpha]*Sz0[alpha];
      }
      for(int j = 0; j < TWOL; j++)
      {
        Szjnow = Szj[j];
        mSzjn = 0;
        for (int alpha = 0; alpha < maxIndexValue; alpha ++)
        {
          mSzjn += eigenvecn[alpha]*eigenvecm[alpha]*Szjnow[alpha];
        }
        //cout << eigenvals[n]-minval << endl;

        for(int b = 0; b < Nb; b++)
          for(int t = 0; t < Nt; t++)
            SzSz[j][b][t] += nSz0m*mSzjn*exponential(-(eigenvals[n]-eigenvals[m])*time[t])*exp(-beta[b]*(eigenvals[n]-minval));
      }
    }
  }
  return SzSz;
}

Matrix<double, Dynamic, Dynamic> Solver::makeSminus(indexstate converttablep)
{
  Matrix<double, Dynamic, Dynamic> Sminusmat(maxIndexValue, TWOL); //Matrix with one row for each basis state in the current mag sector and one column for each site i.

  vector<short int> statevec;
  unsigned int statenum;

  for(int i = 0; i < TWOL; i++)
  {
    for(int alpha = 0; alpha < maxIndexValue; alpha++)
    {
      statevec = statenum_to_statevec(converttable.index_to_state[alpha]);
      if(statevec[i] == +1)
      {
        statevec[i] = -1;
        statenum = statevec_to_statenum(statevec);
        Sminusmat(alpha, i) = converttablep.state_to_index[statenum];
      }
      else Sminusmat(alpha, i) = -1;
    }
  }

  return Sminusmat;
}


Matrix<double, Dynamic, Dynamic> Solver::SminusMATRIX(indexstate converttablep, int siteind)
{
  Matrix<double, Dynamic, Dynamic> Sminusmat(converttablep.MaxIndex, maxIndexValue); //Matrix with one row for each basis state in the current mag sector and one column for each site i.

  vector<short int> statevec;
  unsigned int statenum;

  Sminusmat.setZero();

  for(int alpha = 0; alpha < maxIndexValue; alpha++)
  {
    statevec = statenum_to_statevec(converttable.index_to_state[alpha]);
    if(statevec[siteind] == +1)
    {
      statevec[siteind] = -1;
      statenum = statevec_to_statenum(statevec);
      Sminusmat(converttablep.state_to_index[statenum], alpha) = +1;
    }
  }

  return Sminusmat;
}

vector<vector<vector<complex<double>>>> Solver::SpmCorrMatZERO(Eigen::Matrix<double, -1, 1, 0, -1, 1> eigenvalsp, Matrix<double,Dynamic,Dynamic> eigenvecsp, indexstate converttablep, vector<double> beta, vector<double> time)
{
  Matrix<double,Dynamic,Dynamic> Sminus0eigenbasis;
  vector<Matrix<double,Dynamic,Dynamic>> Sminusjeigenbasis;

  Matrix<double,Dynamic,Dynamic> eigenvecspinv;

  eigenvecspinv = eigenvecsp.inverse();

  for(int j = 0; j < TWOL; j++)
  {
    Sminusjeigenbasis.push_back(eigenvecspinv*SminusMATRIX(converttablep, j)*eigenvecs);
  }
  Sminus0eigenbasis = Sminusjeigenbasis[0];

  double minval = eigenvals[0];
  if(eigenvalsp[0] < minval) minval = eigenvalsp[0];

  vector<vector<complex<double>>> W(Nb, vector<complex<double>>(Nt, zero));
  vector<vector<vector<double>>>partsum(TWOL, vector<vector<double>>(maxIndexValue, vector<double>(converttablep.MaxIndex, 0.0)));
  vector<vector<vector<complex<double>>>>totsum(TWOL, vector<vector<complex<double>>>(Nb, vector<complex<double>>(Nt, zero)));

  int maxInd;
  if(CUTOFF){maxInd = maxIndexValue;}
  else {maxInd = maxIndexValue;}

  for(int n = 0; n < maxInd; n++) //Only for n-sum?
    for(int m = 0; m < converttablep.MaxIndex; m++)
    {
      for(int j = 0; j < TWOL; j++) partsum[j][n][m] += Sminus0eigenbasis(m,n)*Sminusjeigenbasis[j](m,n);

      for(int t = 0; t < Nt; t++)
        for(int b = 0; b < Nb; b++)
        {
          W[b][t] = 0.5*(exp(-beta[b]*(eigenvals[n]-minval))*exponential(-(eigenvals[n]-eigenvalsp[m])*time[t])+exp(-beta[b]*(eigenvalsp[m]-minval))*exponential(-(eigenvalsp[m]-eigenvals[n])*time[t])); //Exponential legger til en i i eksponenten.
          for(int j = 0; j < TWOL; j++) totsum[j][b][t] += W[b][t]*partsum[j][n][m];
        }
      }

  return totsum;
}

vector<vector<vector<complex<double>>>> Solver::SpmCorrMatNN(Eigen::Matrix<double, -1, 1, 0, -1, 1> eigenvalsp, Matrix<double,Dynamic,Dynamic> eigenvecsp, indexstate converttablep, vector<double> beta, vector<double> time)
{
  vector<Matrix<double,Dynamic,Dynamic>> Sminusjeigenbasis;

  Matrix<double,Dynamic,Dynamic> eigenvecspinv;

  eigenvecspinv = eigenvecsp.inverse();

  for(int j = 0; j < TWOL; j++)
  {
    Sminusjeigenbasis.push_back(eigenvecspinv*SminusMATRIX(converttablep, j)*eigenvecs);
  }

  double minval = eigenvals[0];
  if(eigenvalsp[0] < minval) minval = eigenvalsp[0];

  vector<vector<complex<double>>> W(Nb, vector<complex<double>>(Nt, zero));
  vector<vector<vector<double>>>partsum(TWOL, vector<vector<double>>(maxIndexValue, vector<double>(converttablep.MaxIndex, 0.0)));
  vector<vector<vector<complex<double>>>>totsum(TWOL, vector<vector<complex<double>>>(Nb, vector<complex<double>>(Nt, zero)));

  int maxInd;
  if(CUTOFF){maxInd = maxIndexValue;}
  else {maxInd = maxIndexValue;}

  for(int n = 0; n < maxInd; n++) //Only for n-sum?
    for(int m = 0; m < converttablep.MaxIndex; m++)
    {
      for(int j = 0; j < TWOL; j++) partsum[j][n][m] += Sminusjeigenbasis[j](m,n)*Sminusjeigenbasis[(j+1)%TWOL](m,n);

      for(int t = 0; t < Nt; t++)
        for(int b = 0; b < Nb; b++)
        {
          W[b][t] = 0.5*(exp(-beta[b]*(eigenvals[n]-minval))*exponential(-(eigenvals[n]-eigenvalsp[m])*time[t])+exp(-beta[b]*(eigenvalsp[m]-minval))*exponential(-(eigenvalsp[m]-eigenvals[n])*time[t])); //Exponential legger til en i i eksponenten.
          for(int j = 0; j < TWOL; j++) totsum[j][b][t] += W[b][t]*partsum[j][n][m];
        }
      }

  return totsum;
}

vector<vector<vector<complex<double>>>> Solver::SpmCorrMatMID(Eigen::Matrix<double, -1, 1, 0, -1, 1> eigenvalsp, Matrix<double,Dynamic,Dynamic> eigenvecsp, indexstate converttablep, vector<double> beta, vector<double> time)
{
  Matrix<double,Dynamic,Dynamic> SminusMIDeigenbasis;
  vector<Matrix<double,Dynamic,Dynamic>> Sminusjeigenbasis;

  Matrix<double,Dynamic,Dynamic> eigenvecspinv;

  eigenvecspinv = eigenvecsp.inverse();

  for(int j = 0; j < TWOL; j++)
  {
    Sminusjeigenbasis.push_back(eigenvecspinv*SminusMATRIX(converttablep, j)*eigenvecs);
  }
  SminusMIDeigenbasis = Sminusjeigenbasis[int(Nsites/2)];

  double minval = eigenvals[0];
  if(eigenvalsp[0] < minval) minval = eigenvalsp[0];

  vector<vector<complex<double>>> W(Nb, vector<complex<double>>(Nt, zero));
  vector<vector<vector<double>>>partsum(TWOL, vector<vector<double>>(maxIndexValue, vector<double>(converttablep.MaxIndex, 0.0)));
  vector<vector<vector<complex<double>>>>totsum(TWOL, vector<vector<complex<double>>>(Nb, vector<complex<double>>(Nt, zero)));

  int maxInd;
  if(CUTOFF){maxInd = maxIndexValue;}
  else {maxInd = maxIndexValue;}

  for(int n = 0; n < maxInd; n++) //Only for n-sum?
    for(int m = 0; m < converttablep.MaxIndex; m++)
    {
      for(int j = 0; j < TWOL; j++) partsum[j][n][m] += SminusMIDeigenbasis(m,n)*Sminusjeigenbasis[j](m,n);

      for(int t = 0; t < Nt; t++)
        for(int b = 0; b < Nb; b++)
        {
          W[b][t] = 0.5*(exp(-beta[b]*(eigenvals[n]-minval))*exponential(-(eigenvals[n]-eigenvalsp[m])*time[t])+exp(-beta[b]*(eigenvalsp[m]-minval))*exponential(-(eigenvalsp[m]-eigenvals[n])*time[t])); //Exponential legger til en i i eksponenten.
          for(int j = 0; j < TWOL; j++) totsum[j][b][t] += W[b][t]*partsum[j][n][m];
        }
      }

  return totsum;
}

vector<vector<vector<complex<double>>>> Solver::SpmCorr(Matrix<double, Dynamic, Dynamic> Sminusmat, Eigen::Matrix<double, -1, 1, 0, -1, 1> eigenvalsp, Matrix<double,Dynamic,Dynamic> eigenvecsp, indexstate converttablep, vector<double> beta, vector<double> time)
{
  //Compute the correlations in the xy plane: SxiSxj + SyiSyj.

  //Only for one magnetisation sector. All mag sectors need to be added and divided by the partition function.
  //Will scale each term by the minimal energy in the corresponding sector.

  Eigen::Matrix<double, -1, 1, 0, -1, 1> Sminus0 = Sminusmat.col(0);
  //vector<Eigen::Matrix<double, -1, 1, 0, -1, 1>> Sminusj(TWOL);

  vector<vector<complex<double>>> W(Nb, vector<complex<double>>(Nt, zero));
  vector<vector<vector<double>>>partsum(TWOL, vector<vector<double>>(maxIndexValue, vector<double>(converttablep.MaxIndex, 0.0)));
  vector<vector<vector<complex<double>>>>totsum(TWOL, vector<vector<complex<double>>>(Nb, vector<complex<double>>(Nt, zero)));

  double sum0, sumj;

  double minval = eigenvals[0];
  if(eigenvalsp[0] < minval) minval = eigenvalsp[0];

  /*for(int n = 0; n < maxIndexValue; n++)
  {
    for(int m = 0; m < converttablep.MaxIndex; m++)
    {
      //cout << n << "   " << m << endl;

      W = 0.5*(exp(-beta*(eigenvals[n]-minval))*exponential(-(eigenvals[n]-eigenvalsp[m])*t)+exp(-beta*(eigenvalsp[m]-minval))*exponential(-(eigenvalsp[m]-eigenvals[n])*t)); //Exponential legger til en i i eksponenten.

      sumi = 0;
      sumj = 0;

      for(int alpha = 0; alpha < maxIndexValue; alpha++)
      {
        if(Sminusj[alpha] >= 0) sumj += eigenvecs.col(n)[alpha]*eigenvecsp.col(m)[Sminusj[alpha]];
        if(Sminusi[alpha] >= 0) sumi += eigenvecs.col(n)[alpha]*eigenvecsp.col(m)[Sminusi[alpha]];
      }
      totsum += W*sumi*sumj;
    }
  }*/

  for(int n = 0; n < maxIndexValue; n++)
  {
    for(int m = 0; m < converttablep.MaxIndex; m++)
    {
      sum0 = 0;
      for(int alpha = 0; alpha < maxIndexValue; alpha++)
      {
        int beta = Sminus0[alpha];
        if(beta >= 0) sum0 += eigenvecs.col(n)[alpha]*eigenvecsp.col(m)[beta];
      }
      for(int j = 0; j < TWOL; j++)
      {
        sumj = 0;

        for(int alpha = 0; alpha < maxIndexValue; alpha++)
        {
          int beta = Sminusmat.col(j)[alpha];
          if(beta >= 0) sumj += eigenvecs.col(n)[alpha]*eigenvecsp.col(m)[beta];
        }
        partsum[j][n][m] += sum0*sumj;
      }
    }
  }

  for(int n = 0; n < maxIndexValue; n++)
    for(int m = 0; m < converttablep.MaxIndex; m++)
      for(int b = 0; b < Nb; b++)
        for(int t = 0; t < Nt; t++)
        {
          W[b][t] = 0.5*(exp(-beta[b]*(eigenvals[n]-minval))*exponential(-(eigenvals[n]-eigenvalsp[m])*time[t])+exp(-beta[b]*(eigenvalsp[m]-minval))*exponential(-(eigenvalsp[m]-eigenvals[n])*time[t])); //Exponential legger til en i i eksponenten.
          for(int j = 0; j < TWOL; j++) totsum[j][b][t] += W[b][t]*partsum[j][n][m];
        }

  return totsum;
}

double Solver::Sx2(Eigen::Matrix<double, -1, 1, 0, -1, 1> statecoeffs)
{

  //double Solver::Sx2(vector<double> statecoeffs, vector<double> statenumbers)

  //Compute Sx2 for a state with statecoeffs and corresponding statenumbers. I.e. state coeffs gives the coefficients for the different basis states found in statenumbers.

  //I guess I already have the statenumbers if i calculate Sx2 when I am computing the eigenvalues?
  double ans = 0;
  double maxstatenum = (Nh == 0) ? twomax : twomax*(1 + TWOL*TWOLpow[Nh-1]*Nh); //This is an overestimation for Nh>0!
  vector<double> intermediatestates(maxstatenum, 0.0);

  vector<short int> statevec;
  vector<short int> outvec;
  unsigned int outnum;

  for(int i = 0; i < maxIndexValue; i++)
  {
    statevec = statenum_to_statevec(converttable.index_to_state[i]);

    for(int j = 0; j < TWOL; j++)
    {
      outvec = statevec;
      if(statevec[j] == +1 || statevec[j] == -1)
      {
        outvec[j] *= (-1);
        outnum = statevec_to_statenum(outvec);
        intermediatestates[outnum] += statecoeffs[i];
      }
    }
  }

  for(int i = 0; i < intermediatestates.size(); i++)
  {
    ans += intermediatestates[i]*intermediatestates[i];
  }

  return ans*0.25;
}

double Solver::Sy2(Eigen::Matrix<double, -1, 1, 0, -1, 1> statecoeffs)
{
  //Compute Sy2 for a state with statecoeffs and corresponding statenumbers. I.e. state coeffs gives the coefficients for the different basis states found in statenumbers.

  //I guess I already have the statenumbers if i calculate Sy2 when I am computing the eigenvalues?

  double ans = 0;
  double maxstatenum = (Nh == 0) ? twomax : twomax*(1 + TWOL*TWOLpow[Nh-1]*Nh); //This is an overestimation for Nh>0!
  vector<double> intermediatestates(maxstatenum, 0.0);

  vector<short int> statevec;
  vector<short int> outvec;
  unsigned int outnum;

  for(int i = 0; i < maxIndexValue; i++)
  {
    statevec = statenum_to_statevec(converttable.index_to_state[i]);

    for(int j = 0; j < TWOL; j++) //For each site
    {
      outvec = statevec;
      if(statevec[j] == +1 || statevec[j] == -1)
      {
        outvec[j] *= (-1);
        outnum = statevec_to_statenum(outvec);
        intermediatestates[outnum] += -statevec[j]*statecoeffs[i];
      }
    }
  }

  for(int i = 0; i < intermediatestates.size(); i++)
  {
    ans += intermediatestates[i]*intermediatestates[i];
  }

  return ans*0.25;
}

double Solver::Sz2(Eigen::Matrix<double, -1, 1, 0, -1, 1> statecoeffs)
{
  //Compute Sz2 for a state with statecoeffs and corresponding statenumbers. I.e. state coeffs gives the coefficients for the different basis states found in statenumbers.

  //I think this just is summing the squares of the statecoeffs and mulitply by 0.25?
  //How does that make sense. Aren't the statecoeffs normalised?! Yes, it does not make sense, you need to include the total Sz?
  // So for a given nu, Sz2 is just 0.25*(2*nu - Ns)**2

  return 0.25*(2*nu - Ns)*(2*nu - Ns);
}

double Solver::S(Eigen::Matrix<double, -1, 1, 0, -1, 1> statecoeffs)
{
  double ans = Sx2(statecoeffs) + Sy2(statecoeffs) + Sz2(statecoeffs);
  return 0.5*(-1+sqrt(1+4*ans));
}

vector<double> Solver::HoleDens(Eigen::Matrix<double, -1, 1, 0, -1, 1> statecoeffs)
{
  //Compute holedensity for statecoeffs at site i.

  vector<double> ans(TWOL, 0.0);

  unsigned int alphanum;
  double holepos;

  for (int alpha = 0; alpha < maxIndexValue; alpha++)
  {
    //Find the holepositions of state alpha.
    alphanum = converttable.index_to_state[alpha];
    alphanum /= twomax;


    for (int i = 0; i < Nh; i++)
    {
      holepos = alphanum/TWOLpow[Nh-1-i];
      ans[holepos] += statecoeffs[alpha]*statecoeffs[alpha];
      alphanum %= TWOLpow[Nh-1-i];
    }
  }

  return ans;
}


vector<short int> Solver::translate(vector<short int> statevec, int trans)
{
  //Translate statevec by trans to the right
  int veclength = statevec.size();
  vector<short int> transvec(veclength);

  for (int j = 0; j < veclength; j++)
  {
    transvec[(j+trans+veclength)%veclength] = statevec[j];
  }

  return transvec;
}


/*Matrix<double, Dynamic, Dynamic> Solver::TransMat(int trans)
{
  //For each state, translate it by trans to the right and construct the corresponding TransMat.

  Matrix<double, Dynamic, Dynamic> ans(maxIndexValue, maxIndexValue);
  ans.setZero();

  unsigned int statenum, transnum;
  vector<short int> statevec, transvec;
  short unsigned int transind;

  for (int i = 0; i < maxIndexValue; i++)
  {
    statenum = converttable.index_to_state[i];
    statevec = statenum_to_statevec(statenum);
    transvec = translate(statevec, trans);
    transnum = statevec_to_statenum(transvec);
    transind = converttable.state_to_index[transnum];

    ans(transind, i) = 1;
  }

  return ans;
}*/

vector<double> Solver::TransMat(int trans)
{
  //For each state, translate it by trans to the right and construct the corresponding TransMat.

  vector<double> ans(maxIndexValue);

  unsigned int statenum, transnum;
  vector<short int> statevec, transvec;
  short unsigned int transind;

  for (int i = 0; i < maxIndexValue; i++)
  {
    statenum = converttable.index_to_state[i];
    statevec = statenum_to_statevec(statenum);
    transvec = translate(statevec, trans);
    transnum = statevec_to_statenum(transvec);
    transind = converttable.state_to_index[transnum];

    ans[i] = transind;
  }

  return ans;
}


void Solver::GSQvec()
{
  double qstart = clock();
  vector<double> Ta2;

  //Ta1 = TransMat(2); //corresponding to transtaion by a1
  //I don't have to check Ta1. Ta2 i sufficient as Ta1 is (Ta2))^2, so the
  //eigenvalues of Ta1 are the eigenvalues of Ta2 squared.
  Ta2 = TransMat(1); //corresponding to transtaion by a2

  //Now we have the translation matrices. Next, we need to find Q.
  //If the ground state is degenerate, the corresponding eigen states may not be
  //eigensttates of the translation operators, we then need to find the appropriate
  //superpositions. This is done by diagonalising the translation matrices in the
  //ground state subspace.

  //So I should probably first find the GS degeneracy:

  int deg = 1;
  double diff = (eigenvals[deg] - eigenvals[0])/eigenvals[0];

  Qout << "nu = " << nu << "      Lowest energy = " << eigenvals[0] << endl;

  while (abs(diff) < 1e-9 && deg < maxIndexValue)
  {
    deg += 1;
    if (deg == maxIndexValue){break;}
    diff = (eigenvals[deg] - eigenvals[0])/eigenvals[0];
  }
  Qout << "Degeneracy = " << deg << endl;

  //Then construct the subspace matrices:

  //Matrix<double, Dynamic, Dynamic> Ta1subspace(deg, deg);
  Matrix<double, Dynamic, Dynamic> Ta2subspace(deg, deg);
  Ta2subspace.setZero();

  for (int i = 0; i < deg; i++)
    for (int j = 0; j < deg; j++)
    {
      //Ta1subspace(i, j) = eigenvecs.col(i).transpose()*Ta1*eigenvecs.col(j);
      for (int k = 0; k < maxIndexValue; k++) {Ta2subspace(i, j) += eigenvecs.col(i)[Ta2[k]]*eigenvecs.col(j)[k];}
    }

  //Then diagnoalise them:

  //EigenSolver<Matrix<double,Dynamic,Dynamic>> esa2;
  EigenSolver<Matrix<double,Dynamic,Dynamic>> esa2;

  //Eigen::Matrix<complex<double>, -1, 1, 0, -1, 1> eigenvalsa1;
  //Matrix<complex<double>,Dynamic,Dynamic> eigenvecsa1;
  Eigen::Matrix<complex<double>, -1, 1, 0, -1, 1> eigenvalsa2;
  Matrix<complex<double>,Dynamic,Dynamic> eigenvecsa2;

  //esa1 = EigenSolver<Matrix<double,Dynamic,Dynamic>>(Ta1subspace, ComputeEigenvectors);
  esa2 = EigenSolver<Matrix<double,Dynamic,Dynamic>>(Ta2subspace, ComputeEigenvectors);
  //eigenvalsa1 = esa1.eigenvalues();
  //eigenvecsa1 = esa1.eigenvectors();
  eigenvalsa2 = esa2.eigenvalues();
  eigenvecsa2 = esa2.eigenvectors();

  //Then figure out which Q's they correspond to: (loop through all possible Q values?)
  //eigenvalsa1 = e^{iQ.a1} and eigenvalsa2 = e^{iQ.a2}

  //Compute all allowed Q's and all pairs of e^{iQ.a1} and e^{iQ.a2}.
  //For each ground state, find which Q (eigenvalsa1, eigenvalsa2) corresponds to:

  vector<double> possibleQ(2, 0.0);
  vector<vector<double>> Qs(deg, vector<double>(2, 100.0));

  bool found;
  complex<double> diff2;

  for (int i = 0; i < deg; i++)
  {
    found = false;
    for (int ny = 0; ny < 2; ny++)
    {
      for (int nx = 0; nx < Nsites/2; nx++)
      {
        possibleQ = {TWOPI/(Nsites/2.)*nx, TWOPI/sqrt(3)*ny};

        //diff1 = eigenvalsa1[i] - exponential(possibleQ[0]);
        //cout << diff1 << endl;
        diff2 = eigenvalsa2[i] - exponential(0.5*possibleQ[0]+SQRTTHREEOVERTWO*possibleQ[1]);
        //cout << eigenvalsa2[i] << "  " << exponential(0.5*possibleQ[0]+SQRTTHREEOVERTWO*possibleQ[1]) << endl;
        //cout << diff2 << endl;

        if (abs(diff2.real()) < 1e-9 && abs(diff2.imag()) < 1e-9)
        {
          Qs[i] = possibleQ;
          found = true;
          break;
        }
      }
      if (found) break;
    }
    //if (found) cout << "found" << endl;
    //else cout << "not found" << endl;
    //cout << "Q = (" << Qs[i][0] << ", " << Qs[i][1] << ")" << endl;
  }

  //auto inttostringfunc = std::bind(&Solver::myinttostring, this, std::placeholders::_1);

  Eigen::Matrix<complex<double>, -1, 1, 0, -1, 1> Qstate(maxIndexValue);
  for (int i = 0; i < deg; i++)
  {
    Qstate.setZero();
    Qout << "Q = (" << Qs[i][0] << ", " << Qs[i][1] << ")" << endl;
    //for (int k = 0; k < maxIndexValue; k++) Qstate[k] = eigenvecsa2.col(i)[0]*eigenvecs.col(0)[k];
    for (int j = 0; j < deg; j++) for (int k = 0; k < maxIndexValue; k++) {Qstate[k] += eigenvecsa2.col(i)[j]*eigenvecs.col(j)[k];}
    PrintState(Qstate, Qout);
  }

  double qstop = clock();

  Qout << "Q time: " << (qstop-qstart)/CLOCKS_PER_SEC << endl;


  Qout << endl;

  return;
}


vector<short int> Solver::parity(vector<short int> statevec)
{
  //Translate statevec by trans to the right
  int veclength = statevec.size();
  vector<short int> parvec(veclength);

  for (int j = 0; j < veclength; j++)
  {
    parvec[j] = statevec[veclength-1-j];
  }

  return parvec;
}


Matrix<double, Dynamic, Dynamic> Solver::ParMat()
{
  //For each state, translate it by trans to the right and construct the corresponding TransMat.

  Matrix<double, Dynamic, Dynamic> ans(maxIndexValue, maxIndexValue);
  ans.setZero();

  unsigned int statenum, parnum;
  vector<short int> statevec, parvec;
  short unsigned int parind;

  for (int i = 0; i < maxIndexValue; i++)
  {
    statenum = converttable.index_to_state[i];
    statevec = statenum_to_statevec(statenum);
    parvec = parity(statevec);
    parnum = statevec_to_statenum(parvec);
    parind = converttable.state_to_index[parnum];

    ans(parind, i) = 1;
  }

  return ans;
}


void Solver::GSparity()
{
  Matrix<double, Dynamic, Dynamic> P;

  P = ParMat(); //corresponding to transtaion by a2

  //Now we have the parity operator.
  //If the ground state is degenerate, the corresponding eigen states may not be
  //eigensttates of the parity operator, we then need to find the appropriate
  //superpositions. This is done by diagonalising the parity matrix in the
  //ground state subspace.

  //So I should probably first find the GS degeneracy:

  int deg = 1;
  double diff = (eigenvals[deg] - eigenvals[0])/eigenvals[0];

  cout << eigenvals[0] << endl;

  while (abs(diff) < 1e-9 && deg < maxIndexValue)
  {
    deg += 1;
    if (deg == maxIndexValue){break;}
    diff = (eigenvals[deg] - eigenvals[0])/eigenvals[0];
  }
  cout << "Degeneracy = " << deg << endl;

  //Then construct the subspace matrices:

  //Matrix<double, Dynamic, Dynamic> Ta1subspace(deg, deg);
  Matrix<double, Dynamic, Dynamic> Psubspace(deg, deg);

  for (int i = 0; i < deg; i++)
    for (int j = 0; j < deg; j++)
    {
      //Ta1subspace(i, j) = eigenvecs.col(i).transpose()*Ta1*eigenvecs.col(j);
      Psubspace(i, j) = eigenvecs.col(i).transpose()*P*eigenvecs.col(j);
    }

  //Then diagnoalise them:

  //EigenSolver<Matrix<double,Dynamic,Dynamic>> esa2;
  EigenSolver<Matrix<double,Dynamic,Dynamic>> esP;

  //Eigen::Matrix<complex<double>, -1, 1, 0, -1, 1> eigenvalsa1;
  //Matrix<complex<double>,Dynamic,Dynamic> eigenvecsa1;
  Eigen::Matrix<complex<double>, -1, 1, 0, -1, 1> eigenvalsP;
  Matrix<complex<double>,Dynamic,Dynamic> eigenvecsP;

  //esa1 = EigenSolver<Matrix<double,Dynamic,Dynamic>>(Ta1subspace, ComputeEigenvectors);
  esP = EigenSolver<Matrix<double,Dynamic,Dynamic>>(Psubspace, ComputeEigenvectors);
  //eigenvalsa1 = esa1.eigenvalues();
  //eigenvecsa1 = esa1.eigenvectors();
  eigenvalsP = esP.eigenvalues();
  eigenvecsP = esP.eigenvectors();


  //What now? Print parity GS?

  return;
}


/*void Solver::SqSq()
{
  //Compute \ev{S-qSq}

  //Need to choose values of q. For a periodic system it would make sense to only do the allowed q's? How does that work with holes?
  //For OBC, I guess we could choose a fixed set of q's?
  //For OBCx PBCy we could discretize in y-direction?

  int Nx = Nsites/2;
  int Ny = 2;

  double dqx = TWOPI/Nx;
  double dqy = TWOPI/sqrt(3);


  //if (PBC)
  //{
  //  Nx = Nsites/2;
  //}
  //else
  //{
  //  Nx = 20;
  //}

  int Nq = Nx*Ny;

  vector<double> q(2);

  for (int nqy = 0; nqy < 2; nqy++)
    for (int nqx = 0; nqx < Nx; nqx++)
    {
      q = {dqx*nqx, dqy*nqy};

      //THIS IS WHERE YOU ARE WORKING CURRENTLY

      //Compute SqSq and write to file q SqSq. Cannot compute it here, it must be computed for each sector and added properly together.
    }

  return;
}


void Solver::LocMag();
{
  //Compute the local magnetisation.
  //Maybe the expression for the local magnetisation should depend on some variable?
  //Start with the AFM on triangular lattice and extend to other cases later.

  //Expectation value of local magnetisation?

  //Must compute for each sector separately.
}*/


string Solver::index_to_string(unsigned short int stateind)
{
  unsigned int statenum = converttable.index_to_state[stateind];
  vector<short int> statevec = statenum_to_statevec(statenum);

  string statestring;

  if (statevec[0] == -1) statestring = "-";
  else if (statevec[0] == 0) statestring = "o";
  else if (statevec[0] == +1) statestring = "+";

  for (int pos = 1; pos < Nsites; pos++)
  {
    if (statevec[pos] == -1) statestring += "-";
    else if (statevec[pos] == 0) statestring += "o";
    else if (statevec[pos] == +1) statestring += "+";
  }

  return statestring;
}


void Solver::PrintState(Eigen::Matrix<complex<double>, -1, 1, 0, -1, 1> state, ostream& os)
{
  const int WORDLENGTH = index_to_string(0).size();
  const int NWORDSONALINE = PAGEWIDTH/(WORDLENGTH+SPACING+PRESPACING);

  PrettyPrint pp(os,NWORDSONALINE,WORDLENGTH,SPACING);

  vector<cspair> estate(maxIndexValue);
  for(int i = 0; i < maxIndexValue; i++)
  {
    estate[i] = cspair(state[i],i);
  }

  sort(estate.begin(),estate.end(),bigabsfirst);

  double startcoeff = abs(estate[0].c);
  double currentcoeff = startcoeff;
  int i = 0;

  while(i < maxIndexValue && abs(estate[i].c) > COEFFLIMIT*startcoeff)
  {
    os << "+" << currentcoeff << "(";
    while(abs(abs(estate[i].c) - currentcoeff) < 1e-10)
    {
      pp.add(estate[i].c, index_to_string(estate[i].s));
	    i++;
	  }
	  pp.flush();
	  if(i < maxIndexValue) {currentcoeff = abs(estate[i].c);}
	  os << ")" << endl;
  }
}


void Solver::WriteEigvals()
{
  ofstream Outfile(dir + "eigvals.txt", std::ios_base::app);
  if (!Outfile.is_open())
     cout<<"Could not open file" << endl;

  Outfile.precision(17);
  Outfile << nu << "      " << eigenvals.transpose() << "\n";
}


void Solver::WriteSzStot()
{
  ofstream Outfile(dir + "SzStot.txt", std::ios_base::app);
  if (!Outfile.is_open())
     cout<<"Could not open file" << endl;

  Outfile.precision(17);
  for(int i = 0; i < maxIndexValue; i++)
  {
    Outfile << 0.5*(2*nu - Ns) << "   " << S(eigenvecs.col(i)) << endl;
  }
  Outfile << "\n";
}

void Solver::WriteCorr(vector<double> beta, vector<double> time, vector<vector<vector<complex<double>>>> z, vector<vector<vector<complex<double>>>> pm, string file)
{
  cout << "Starter her" << endl;

  ofstream Outfile(dir + file, std::ios_base::app);
  if (!Outfile.is_open())
     cout<<"Could not open file" << endl;

  cout << "Inni her" << endl;

  Outfile.precision(17);
  for(int b = 0; b < beta.size(); b++)
    for(int t = 0; t < time.size(); t++)
    {
      Outfile << beta[b] << "   " << time[t] << "   ";
      for(int i = 0; i < z.size(); i++) Outfile << z[i][b][t] << "   ";
      for(int i = 0; i < pm.size(); i++) Outfile << pm[i][b][t] << "   ";
      Outfile << endl;
    }
}

void Solver::WritePartition(vector<double> beta, vector<double> partition)
{
  ofstream Outfile(dir + "Partition.txt", std::ios_base::app);
  if (!Outfile.is_open())
     cout<<"Could not open file" << endl;

  for(int b = 0; b < beta.size(); b++)
  {
    Outfile << beta[b] << "   " << partition[b] << endl;
  }
}

void Solver::WriteHoleDens()
{
  //Find and write to file the expectation value of the holedensity at different sites for the ground state of each sector?

  //Need a function which computes the expectation value.

  vector<double> HoleDensvec = HoleDens(eigenvecs.col(0));

  //Then do the writing to file here.

  ofstream Outfile(dir + "GSHoleDensity.txt", std::ios_base::app);
  if (!Outfile.is_open())
     cout<<"Could not open file" << endl;

  Outfile << nu << "   " << eigenvals[0] << "   ";

  for (int i = 0; i < TWOL; i++)
  {
    Outfile << HoleDensvec[i] << "   ";
  }

  Outfile << endl;
}

void Solver::WriteHoleCorr(vector<double> beta, vector<double> time, vector<vector<vector<complex<double>>>> NhNh, string file)
{

  ofstream Outfile(dir + file, std::ios_base::app);
  if (!Outfile.is_open())
     cout<<"Could not open file" << endl;

  Outfile.precision(17);
  for(int b = 0; b < beta.size(); b++)
    for(int t = 0; t < time.size(); t++)
    {
      Outfile << beta[b] << "   " << time[t] << "   ";
      for(int i = 0; i < NhNh.size(); i++) Outfile << NhNh[i][b][t] << "   ";
      Outfile << endl;
    }
}


void Solver::resetdatafiles()
{
  ofstream EigFile(dir + "eigvals.txt");
  if (!EigFile.is_open())
     cout<<"Could not open file" << endl;

  ofstream SFile(dir + "SzStot.txt");
  if (!SFile.is_open())
    cout<<"Could not open file" << endl;

  ofstream PartitionFile(dir + "Partition.txt");
  if (!PartitionFile.is_open())
    cout<<"Could not open file" << endl;

  ofstream CorrFile(dir + "Corr.txt");
  if (!CorrFile.is_open())
    cout<<"Could not open file" << endl;

  ofstream CorrFileNN(dir + "CorrNN.txt");
  if (!CorrFileNN.is_open())
    cout<<"Could not open file" << endl;

  ofstream CorrFileMID(dir + "CorrMID.txt");
  if (!CorrFileMID.is_open())
    cout<<"Could not open file" << endl;

  ofstream GSHoleDensFile(dir + "GSHoleDensity.txt");
  if (!GSHoleDensFile.is_open())
    cout<<"Could not open file" << endl;

  ofstream HoleCorr(dir + "HoleCorr.txt");
  if (!HoleCorr.is_open())
    cout<<"Could not open file" << endl;

  ofstream HoleCorrNN(dir + "HoleCorrNN.txt");
  if (!HoleCorrNN.is_open())
    cout<<"Could not open file" << endl;

  ofstream HoleCorrMID(dir + "HoleCorrMID.txt");
  if (!HoleCorrMID.is_open())
    cout<<"Could not open file" << endl;
}

#endif
