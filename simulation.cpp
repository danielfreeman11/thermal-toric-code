#include <vector>
#include "lattice.h"
#include <math.h>
#include <stdlib.h>

using namespace std;

using std::vector;
using std::tuple;


float AdvanceCounter(float c, float curt, float maxt, vector<float>* DiscObs, vector<tuple<float,int,int,int>>* Obs)
{
	if (curt > (c / 10.) * maxt)
	{
		//cout << c << "\t" << curt << "\t" << maxt << "\t" << (c / 10.) * maxt << endl;
		if ((get<1>(Obs->at(Obs->size()-2)) != 0 || get<2>(Obs->at(Obs->size()-2)) != 0) && c <=10 ) //Pick out the second-to last observed winding operators
		{																							 //This is because this loop checks for the first time that the simulation time has exceeded the bin.  We also have to account for a single timestep being larger than a single bin.  Thus, function checks for the first time that exceeds the bin number, adds the observable at the previous time to the bin (because it must be the observable that is observed at time t = c/10. * maxtime), and then advances by one bin.  On the next pass, it does the same thing until the bin number has advanced enough.
			//cout << "Writing observable: " << get<1>(Obs->at(Obs->size()-2)) << "\t" << get<2>(Obs->at(Obs->size()-2)) << " at time " << get<0>(Obs->at(Obs->size()-2)) << endl;
			(DiscObs)->at((int) c) = (DiscObs)->at((int) c) + 1.;
		}
		
		
		return 1 + AdvanceCounter(c + 1, curt, maxt, DiscObs, Obs);
		
		//write Obs at immediately preceding index to DiscreteObs
		//advance c by 1 + AdvanceCounter
	}
	else
	{
		return 0.;
	}
	

}



int main ()
{

	srand (time(NULL));

	vector<float> Temp = {.05, .06, .07, .08};
	vector<float> MeanLifetimes;




	//first, need a decent estimate of the lifetime
	for (int jj = 0; jj < 1; jj++)
	{
		float cc = exp(-1./Temp[jj]);
		vector<float> *Lifetimes = new vector<float>;
		for (int kk = 0; kk < 10; kk++)
		{

			Lattice *myLattice = new Lattice(32,cc,1.0,Temp[jj]);
			/*if (kk==0)
			{
				myLattice->Print_Lattice();
				//cout << Temp[jj];
			}*/
			
			vector<tuple<float,int,int,int>> *Observables = new vector<tuple<float,int,int,int>>;
			tuple<float,int,int,int> o1 (0.,0,0,0);
			Observables->push_back(o1);

			while ((myLattice->SectorX == 0 && myLattice->SectorY == 0) || myLattice->anyons.size() > 0)
			{	
				
				myLattice->Update();
				myLattice->Update_Sector();
				//myLattice->Print_Lattice();
				
				tuple<float,int,int,int> o(myLattice->TotalTime,myLattice->SectorX,myLattice->SectorY,myLattice->anyons.size());
				//cout << get<0>(o) << "\t" << get<1>(o) << "\t" <<get<2>(o) << "\t" << myLattice->anyons.size() << endl;
				/*if (myLattice->anyons.size() == 2)
				{
					cout << sqrt(pow((myLattice->anyons[0]->x-myLattice->anyons[1]->x),2) + pow((myLattice->anyons[0]->y-myLattice->anyons[1]->y),2)) << endl;
				}*/
				Observables->push_back(o);
			}

			Lifetimes->push_back(get<0>(Observables->back()));
			
			delete Observables;
			delete myLattice;
			
			//cout << Lifetimes->back() << endl;
			

		}
		double sum = accumulate(Lifetimes->begin(), Lifetimes->end(), 0.0);
		double mean = sum / Lifetimes->size();

		double sq_sum = inner_product(Lifetimes->begin(), Lifetimes->end(), Lifetimes->begin(), 0.0);
		double stdev = sqrt(sq_sum / Lifetimes->size() - mean * mean);
		
		delete Lifetimes;

		cout << Temp[jj] << "\t" << mean << "\t" << stdev << endl;
		MeanLifetimes.push_back(mean);
		
	}







	for (int jj = 0; jj < 1; jj++)
	{
		float cc = exp(-1./Temp[jj]);
		vector<float> DiscreteObs = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
		for (int kk = 0; kk < 1000; kk++)
		{

			Lattice *myLattice = new Lattice(32,cc,1.0,Temp[jj]);
			/*if (kk==0)
			{
				myLattice->Print_Lattice();
				//cout << Temp[jj];
			}*/
			
			vector<tuple<float,int,int,int>> *Observables = new vector<tuple<float,int,int,int>>;
			tuple<float,int,int,int> o1 (0.,0,0,0);
			Observables->push_back(o1);
			
			int counter = 1; //Keeps track of which bin of DiscreteObs to write into

			while (myLattice->TotalTime < MeanLifetimes[jj])
			{	
				
				myLattice->Update();
				myLattice->Update_Sector();
				//myLattice->Print_Lattice();
				
				tuple<float,int,int,int> o(myLattice->TotalTime,myLattice->SectorX,myLattice->SectorY,myLattice->anyons.size());
				//cout << get<0>(o) << "\t" << get<1>(o) << "\t" <<get<2>(o) << "\t" << myLattice->anyons.size() << endl;
				/*if (myLattice->anyons.size() == 2)
				{
					cout << sqrt(pow((myLattice->anyons[0]->x-myLattice->anyons[1]->x),2) + pow((myLattice->anyons[0]->y-myLattice->anyons[1]->y),2)) << endl;
				}*/
				Observables->push_back(o);
				
				counter = counter + (int)(AdvanceCounter((float) counter, myLattice->TotalTime, MeanLifetimes[jj], &DiscreteObs, Observables));

				
				
			}

			
			delete Observables;
			delete myLattice;
			
			//cout << Lifetimes->back() << endl;
			

		}
		cout << "Observable histogram: " << endl;
		for (int a = 0; a < 10; a++)
			cout << (((float)a) * MeanLifetimes[jj]/10.) << "\t" << DiscreteObs[a] << endl;
		
	}








	//Anyon* aa = new Anyon(4,5);
	//Anyon* b = new Anyon(4,4);

	//myLattice->Add_Anyon(aa);

	//myLattice->Print_Lattice();

	//myLattice->Add_Anyon(b);

	//myLattice->Print_Lattice();*/

	/*for (int i = 0; i < Observables.size(); i++)
	{

		cout << get<0>(Observables[i]) << "\t" << get<1>(Observables[i]) << "\t" << get<2>(Observables[i]) << "\t" << get<3>(Observables[i]) <<  endl;

	}*/
	/*myLattice->Update();
	myLattice->Print_Lattice();
	myLattice->Update();
	myLattice->Print_Lattice();
	myLattice->Update();
	myLattice->Print_Lattice();
	myLattice->Update();
	myLattice->Print_Lattice();
	myLattice->Update();

	myLattice->Print_Lattice();*/


	return 0;
}