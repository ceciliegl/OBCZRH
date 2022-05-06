#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<algorithm>

#include <complex>

#include <Eigen/Dense>
using namespace Eigen;

using namespace std;


const int PAGEWIDTH=120;
const int SPACING=0;    // extra spaces between words
const int PRESPACING=5; // spaces needed for printing the sign (use 4)
const double COEFFLIMIT=0.5; // only print out those with coeffs > COEFFLIMIT times the biggest one.

class PrettyPrint
{
public:
  PrettyPrint(ostream& os_in,int in_W,int in_ssize,int in_spacing):os(os_in),W(in_W),ssize(in_ssize),spacing(in_spacing),prespacing(PRESPACING),counter(0),cof(W),str(W){}
  ~PrettyPrint(){flush();}
  void add(complex<double> a,string s)
  {
    cof[counter]=a; str[counter]=s;
    counter++;
    if(counter==W){flush();}
  }
  void flush()  // causes a printout, this can be also be called to force printout even if counter<W.
  {
    print();

    //reset:
    for(int i=0; i<W; i++){cof[i]=0.; str[i]="";}
    counter=0;
  }

  void print()     // the actual print routine
  {
    std::ios oldState(nullptr);
    oldState.copyfmt(std::cout);

    if(counter==0) return;

    os << endl;
    for(int w=0; w<counter; w++)
    {
      // odd indices
      os << std::string(prespacing,' ');
      for(int i=1; i<ssize; i+=2){os << " " << str[w][i];}
      os << std::string(spacing,' ');
    }
    os << endl;

    for(int w=0; w<counter; w++)
    {
      //print the sign of the state:
      //os << (cof[w]<0 ? " m  ": " p  ");
      double angle = arg(cof[w]);
      os << " " << fixed << setprecision(2) << (angle < 0 ? angle + TWOPI : angle); // the length of these strings must equal prespacing
      // even indices
      for(int i=0; i<ssize; i+=2){os << str[w][i] << " ";}
      os << std::string(spacing,' ');
    }
    os << endl;

    std::cout.copyfmt(oldState);
  }

private:
  ostream& os;
  const int W;         // number of words/states on a line, before flushing
  const int ssize;      // length of each string state.
  const int spacing;    // separation between states on a line.
  const int prespacing; // for use with printing the sign
  int counter;
  vector<complex<double>> cof; // container for coefficients
  vector<string> str; // container for string defining states
};

struct cspair
{
public:
  cspair(complex<double> c_in=0,int s_in=0):c(c_in),s(s_in){}
  complex<double> c;
  int s;
};

bool bigabsfirst(cspair a,  cspair b)
{
  if(abs(abs(a.c)-abs(b.c)) < 1e-13){return a.s<b.s;}
  return abs(a.c)>abs(b.c);
}


/*
class MakePretty
{
  friend ostream& operator<<(ostream& os,const MakePretty& mp){ mp.Print(os); return os; }
 public:
 MakePretty(Eigen::DenseBase<Eigen::Matrix<double,-1,-1,0,-1,-1> >::ColXpr& in_state,string (*func)(int)):state(in_state),inttostring(func),WORDLENGTH(inttostring(0).size()),NWORDSONALINE(PAGEWIDTH/(WORDLENGTH+SPACING+PRESPACING)){}
  void const Print(ostream& os) const
  {
    PrettyPrint pp(os,NWORDSONALINE,WORDLENGTH,SPACING);

    int N=state.size();
    vector<cspair> estate(N);
    for(int i=0; i<N; i++)
      {
	estate[i]=cspair(state[i],i);
      }
    sort(estate.begin(),estate.end(),bigabsfirst);

    double startcoeff=abs(estate[0].c);
    double currentcoeff=startcoeff;
    int i=0;
    while(i<N && abs(estate[i].c) > COEFFLIMIT*startcoeff)
      {
	os << "+" << currentcoeff << "(";
	while( abs(estate[i].c)==currentcoeff)
	  {
	    pp.add(estate[i].c,inttostring(estate[i].s));
	    i++;
	  }
	pp.flush();
	if(i<N){ currentcoeff=abs(estate[i].c);}
	os << ")" << endl;
      }
  }
 private:
  Eigen::DenseBase<Eigen::Matrix<double,-1,-1,0,-1,-1> >::ColXpr& state;
  string (*inttostring)(int);
  const int WORDLENGTH;
  const int NWORDSONALINE;
};
*/
