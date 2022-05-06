#include "makepretty.h"

string myinttostring(int a)
{
  return string("-o+-++-+");
}



int main()
{
  // make a test vector
  int N=40;
  vector<double> state(N);
  
  for(int i=0; i<N; i++)
    {
      state[i]=0.1*(rand()%5)-0.2;
    }

  //usage:
  cout << MakePretty(state,myinttostring);
}



