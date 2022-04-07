#include <fstream>
#include<string>

using namespace std;

class ReadInputFiles
{
public:
  string filename_parameters;

  int Nsites, Nh;
  double tl, tr, Jzl, Jzr, Jpml, Jpmr;

  bool OBC;

  bool PYROCHLORE;

  bool EIGVECS, ZEROCORR, NNCORR, MIDCORR;
  bool RESETFILES;

  ReadInputFiles(string filename_param);
  void generate();
};

ReadInputFiles::ReadInputFiles(string filename_param)
{
  filename_parameters = filename_param;
}



void ReadInputFiles::generate()
{
  ifstream parameters(filename_parameters);
  if (!parameters)
  {
    cerr << "Unable to open file " << filename_parameters << "." << endl;
    exit(1);   // call system to stop
  }

  string nothing;
  parameters >> nothing >> nothing >> Nsites;
  parameters >> nothing >> nothing >> Nh;
  parameters >> nothing >> nothing >> tl; parameters >> nothing >> nothing >> tr; parameters >> nothing >> nothing >> Jzl; parameters >> nothing >> nothing >> Jzr; parameters >> nothing >> nothing >> Jpml; parameters >> nothing >> nothing >> Jpmr;
  parameters >> nothing >> nothing >> OBC;
  parameters >> nothing >> nothing >> PYROCHLORE;
  parameters >> nothing >> nothing >> EIGVECS;
  parameters >> nothing >> nothing >> ZEROCORR;
  parameters >> nothing >> nothing >> NNCORR;
  parameters >> nothing >> nothing >> MIDCORR;
  parameters >> nothing >> nothing >> RESETFILES;

  parameters.close();
}
