#include <vector>
#include "lattice.h"
#include <math.h>
#include <stdlib.h>

using namespace std;

using std::vector;
using std::tuple;


float sim_fidelity = 1000.;


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

float AdvanceCounterEntanglementFine(float c, float curt, float maxt, vector<float>* DiscObs, vector<tuple<float,int,int,int,int**,int**>>* Obs)
{
	if (curt > (c / sim_fidelity) * maxt)
	{
		//cout << c << "\t" << curt << "\t" << maxt << "\t" << (c / 10.) * maxt << endl;
		if ((get<1>(Obs->at(Obs->size()-2)) != 0 || get<2>(Obs->at(Obs->size()-2)) != 0) && c <=sim_fidelity ) //Pick out the second-to last observed winding operators
		{																							 //This is because this loop checks for the first time that the simulation time has exceeded the bin.  We also have to account for a single timestep being larger than a single bin.  Thus, function checks for the first time that exceeds the bin number, adds the observable at the previous time to the bin (because it must be the observable that is observed at time t = c/10. * maxtime), and then advances by one bin.  On the next pass, it does the same thing until the bin number has advanced enough.
			//cout << "Writing observable: " << get<1>(Obs->at(Obs->size()-2)) << "\t" << get<2>(Obs->at(Obs->size()-2)) << " at time " << get<0>(Obs->at(Obs->size()-2)) << endl;
			(DiscObs)->at((int) c) = (DiscObs)->at((int) c) + 1.;
		}
		
		
		return 1 + AdvanceCounterEntanglementFine(c + 1, curt, maxt, DiscObs, Obs);
		
		//write Obs at immediately preceding index to DiscreteObs
		//advance c by 1 + AdvanceCounter
	}
	else
	{
		return 0.;
	}



}

/*float AdvancePair(float c1, float c2, float curt1, float curt2, vector<float>* DiscObs, vector<tuple<float,int,int,int,int**,int**>>* Obs1, vector<tuple<float,int,int,int,int**,int**>>* Obs2)
{
	//Updates DiscObs with measurements occurring at the same time for observables obs1 and obs2
	if (curt > (c / sim_fidelity) * maxt)
	{
		//cout << c << "\t" << curt << "\t" << maxt << "\t" << (c / 10.) * maxt << endl;
		if ((get<1>(Obs->at(Obs->size()-2)) != 0 || get<2>(Obs->at(Obs->size()-2)) != 0) && c <=sim_fidelity ) //Pick out the second-to last observed winding operators
		{																							 //This is because this loop checks for the first time that the simulation time has exceeded the bin.  We also have to account for a single timestep being larger than a single bin.  Thus, function checks for the first time that exceeds the bin number, adds the observable at the previous time to the bin (because it must be the observable that is observed at time t = c/10. * maxtime), and then advances by one bin.  On the next pass, it does the same thing until the bin number has advanced enough.
			//cout << "Writing observable: " << get<1>(Obs->at(Obs->size()-2)) << "\t" << get<2>(Obs->at(Obs->size()-2)) << " at time " << get<0>(Obs->at(Obs->size()-2)) << endl;
			(DiscObs)->at((int) c) = (DiscObs)->at((int) c) + 1.;
		}
		
		
		return 1 + AdvanceCounterEntanglementFine(c + 1, curt, maxt, DiscObs, Obs);
		
		//write Obs at immediately preceding index to DiscreteObs
		//advance c by 1 + AdvanceCounter
	}
	else
	{
		return 0.;
	}
	





}*/


void ProcessObservables(int Current_Index, vector<float>* DiscObs, tuple<float,int,int,int,int**,int**> Obs1, tuple<float,int,int,int,int**,int**> Obs2)
{

	if (get<4>(Obs1)[0][0]==1 && get<4>(Obs2)[0][0]==1)
	{
		(DiscObs)->at(Current_Index) = (DiscObs)->at(Current_Index) + 1.;
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
			cout << Observables->size() << endl;
			
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
		//vector<float> DiscreteObs((int)sim_fidelity,0.);// = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
		
		//vector<float> DiscreteObs = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
		std::vector<float> DiscreteObs( ((int) sim_fidelity)+1, 0.0);
		for (int kk = 0; kk < 100; kk++)
		{

			Lattice *myLattice1 = new Lattice(32,cc,1.0,Temp[jj]);
			Lattice *myLattice2 = new Lattice(32,cc,1.0,Temp[jj]);
			/*if (kk==0)
			{
				myLattice->Print_Lattice();
				//cout << Temp[jj];
			}*/
			
			vector<tuple<float,int,int,int,int**,int**>> *Observables1 = new vector<tuple<float,int,int,int,int**,int**>>;
			vector<tuple<float,int,int,int,int**,int**>> *Observables2 = new vector<tuple<float,int,int,int,int**,int**>>;

			tuple<float,int,int,int, int**,int**> o1 (0.,0,0,0,myLattice1->error_map,myLattice1->lattice_verts);
			tuple<float,int,int,int, int**,int**> o2 (0.,0,0,0,myLattice2->error_map,myLattice2->lattice_verts);
			Observables1->push_back(o1);
			Observables2->push_back(o2);
			
			bool Past_Discrete_Observation_Time_1 = false;
			bool Past_Discrete_Observation_Time_2 = false;
			
			int Observables_Index_1 = 0;
			int Observables_Index_2 = 0;
			
			int Discrete_Observation_Index = 1; //Keeps track of which bin of DiscreteObs to write into
			
			
			while ((myLattice1->TotalTime < MeanLifetimes[jj]) && (myLattice1->TotalTime < MeanLifetimes[jj]))
			{
				//Alternates simulation so that they're close to each other in time.  This allows calculation of observables to be easier.
				if (myLattice1->TotalTime <= myLattice2->TotalTime)
				{
					myLattice1->Update();
					myLattice1->Update_Sector();
					tuple<float,int,int,int, int**,int**> o(myLattice1->TotalTime,myLattice1->SectorX,myLattice1->SectorY,myLattice1->anyons.size(),myLattice1->error_map, myLattice1->lattice_verts);
					Observables1->push_back(o);
				}
				else
				{
					myLattice2->Update();
					myLattice2->Update_Sector();
					tuple<float,int,int,int, int**,int**> o(myLattice2->TotalTime,myLattice2->SectorX,myLattice2->SectorY,myLattice2->anyons.size(),myLattice2->error_map, myLattice2->lattice_verts);
					Observables2->push_back(o);
				}
				
				if ((myLattice1->TotalTime > ((float)Discrete_Observation_Index / sim_fidelity) * MeanLifetimes[jj]) && Past_Discrete_Observation_Time_1 == false) //Checks to see if the first lattice has exceeded the current discrete observation window for the first time.  If so, sets the boolean to true and stores the index of the relevant observation
				{
					Past_Discrete_Observation_Time_1 = true;
					Observables_Index_1 = ((Observables1->size())-2);
				}
				
				if ((myLattice2->TotalTime > ((float)Discrete_Observation_Index / sim_fidelity) * MeanLifetimes[jj]) && Past_Discrete_Observation_Time_2 == false) //Checks to see if the second lattice has exceeded the current discrete observation window for the first time.  If so, sets the boolean to true and stores the index of the relevant observation
				{
					Past_Discrete_Observation_Time_2 = true;
					Observables_Index_2 = ((Observables2->size())-2);
				}
				
				if (Past_Discrete_Observation_Time_1 == true && Past_Discrete_Observation_Time_2 == true)
				{
					//Make joint measurement
					//But also make sure that you don't skip over an observation window.
					
					while ( (myLattice1->TotalTime > (Discrete_Observation_Index / sim_fidelity) * MeanLifetimes[jj]) && (myLattice1->TotalTime > (Discrete_Observation_Index / sim_fidelity) * MeanLifetimes[jj]) && Discrete_Observation_Index < sim_fidelity) //This while loop makes sure that the discrete observables list is filled in correctly.  It's possible that an update could occur that's larger than a single discrete_observation window for *both* simulations.  Thus, the discrete observation would be the same for two samples in a row (because nothing changed in the simulation for two discrete observation cycles).  Once the discrete observation window lands after one of either of the total simulation times (or after both), the discrete observation index stops being iterated.
					{	
						
						//cout << "Bad things happening" << Discrete_Observation_Index << "\t" << Observables_Index_1 << "\t" << Observables_Index_2 << "\t" <<  endl;
						//cout << "Bad things happening1" << &Observables1->at(Observables_Index_1) <<  endl;
						//cout << "Bad things happening1" << &Observables2->at(Observables_Index_2) <<  endl;
						ProcessObservables(Discrete_Observation_Index, &DiscreteObs, Observables1->at(Observables_Index_1), Observables2->at(Observables_Index_2)); //Do the thing
						
						//cout << "Bad things happening" << Discrete_Observation_Index <<  endl;
						
						Discrete_Observation_Index += 1;
						
						//Checks to see if the next discrete observation window is before or after the current simulation time.  Updates the boolean accordingly.
						if (myLattice1->TotalTime < ((float)Discrete_Observation_Index / sim_fidelity) * MeanLifetimes[jj])
						{
							Past_Discrete_Observation_Time_1 = false;
						}
						
						if (myLattice2->TotalTime < ((float)Discrete_Observation_Index / sim_fidelity) * MeanLifetimes[jj])
						{
							Past_Discrete_Observation_Time_2 = false;
						}
						
						
						
						
					}
				}
				
				//cout << "Core dump2?" << Discrete_Observation_Index << endl;
				
				/*
				if mylattice1time > discobstime for first time since iterating:
					store this index_1
					
				if mylattice2time > discobstime for first time since iterating
					store this index_2
					
				if both index_1 and index_2 have been advanced past the discreteobstime
					store joint observable at index_1 - 1 and index_2 -1
					advance discobstime
					
					reset first time flag for both.*/
			
			}
			
			cout << "Core dump?" << kk << endl;

			/*while (myLattice1->TotalTime < MeanLifetimes[jj])
			{	
				
				myLattice1->Update();
				myLattice1->Update_Sector();		
				tuple<float,int,int,int, int**,int**> o(myLattice1->TotalTime,myLattice1->SectorX,myLattice1->SectorY,myLattice1->anyons.size(),myLattice1->error_map, myLattice1->lattice_verts);
				Observables1->push_back(o);
				counter = counter + (int)(AdvanceCounterEntanglementFine((float) counter, myLattice1->TotalTime, MeanLifetimes[jj], &DiscreteObs, Observables1));
			}
			
			while (myLattice2->TotalTime < MeanLifetimes[jj])
			{	
				myLattice2->Update();
				myLattice2->Update_Sector();
				tuple<float,int,int,int, int**,int**> o(myLattice2->TotalTime,myLattice2->SectorX,myLattice2->SectorY,myLattice2->anyons.size(),myLattice2->error_map, myLattice2->lattice_verts);
				Observables2->push_back(o);
				//counter = counter + (int)(AdvanceCounterEntanglementFine((float) counter, myLattice1->TotalTime, MeanLifetimes[jj], &DiscreteObs, Observables1));
			}*/ 

			
			delete Observables1;
			delete Observables2;
			
			delete myLattice1;
			delete myLattice2;
			
			//cout << Lifetimes->back() << endl;
			

		}
		cout << "Observable histogram: " << endl;
		for (int a = 0; a < sim_fidelity; a++)
			if (a % int((sim_fidelity / 10.)) == 0)
			{
				cout << (((float)a) * MeanLifetimes[jj]/10.) << "\t" << DiscreteObs[a] << endl;
			}
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
