#include <list>
#include <iostream>
#include <algorithm>
#include <math.h>

using namespace std;

int mod(int a, int b)
{ return (a%b+b)%b; }

class Anyon
{
	public:
		int x;
		int y;
		Anyon(int, int);

};

Anyon::Anyon(int xx, int yy)
{
	x = xx;
	y = yy;
}


class Edge
{

	public:
		Edge() {};
		Edge(int, int, int);
		int edgetype; //0=creatable, 1=translateable, 2=annihilateable
		int x;
		int y12; //weird name to reemphasize that there are twice as many y indices as x indices because of how the counting on the lattice happes
		void Change_Type(int);
		void Update_Anyons(Anyon* a, int m); //m = 0 for remove, m = 1 for add
		Anyon *A1 = NULL; //If swapable, only A1 has an anyon referenced
		Anyon *A2 = NULL; //Otherwise, if killable both A1 and A2 are referenced
		

};


Edge::Edge(int etype, int xx, int yy)
{
	edgetype = etype;
	x = xx;
	y12 = yy;
	
}

void Edge::Change_Type(int etype)
{
	edgetype = etype;
}

void Edge::Update_Anyons(Anyon * a, int m)
{
	if (m==0)
	{
		if (A1 && A2) //if both anyons are already occupied, we delete the matching case and reshuffled A2 to A1
		{
			if (a == A1)
			{
				A1 = A2;
				A2 = NULL;
			}
			else
			{
				A2 = NULL; //else A2 must be 'a' thus it must be deleted
			}
		}
		else //else there's either one anyon or no anyons.  in either case, it's safe to null it because it must be 'a'.
		{
			A1=NULL;
		}
		
	}
	else 
	{
		if (A1) //When adding anyons, either the first one is taken already, or it isn't
		{
			A2 = a;
		}
		else
		{
			A1 = a;
		}
	}
	

	

}

class Lattice
{
	public:
		int ArraySize;
		int SectorX = 0;
		int SectorY = 0;
		
		float TotalTime = 0.;
		
		float cRate;
		float kRate;
		float sRate;
		
		int **lattice_verts;
		int **error_map;
		Edge **lattice_edges;
		vector<Anyon *> anyons;
		vector<Edge *> creatable;
		vector<Edge *> killable;
		vector<Edge *> swapable;
		
		Lattice(int,float,float,float);
		~Lattice();
		
		void Add_Anyon(Anyon* a, bool u);
		void Remove_Anyon(Anyon* a, bool u);
		void Translate_Anyon(Anyon* a, int newx, int newy);
		void Print_Lattice();
		void Update();
		void Update_Sector();

};
	
Lattice::Lattice(int size, float c, float k, float s)
	{
	ArraySize = size;
	cRate = c;
	kRate = k;
	sRate = s;
	
	lattice_verts = new int*[size];
	for (int i =0; i < size; i++)
	{
		lattice_verts[i] = new int[size];
	}
	
	lattice_edges = new Edge*[size];
	for (int i =0; i < size; i++)
	{
		lattice_edges[i] = new Edge[2*size];
	}
	
	error_map = new int*[size];
	for (int i =0; i < size; i++)
	{
		error_map[i] = new int[2*size];
	}
	
	int l = size;
	
	for (int i=0; i<l; i++)
		for (int j=0; j<l; j++)
			{
			lattice_verts[i][j]=0;
			}
	
	for (int i=0; i<l; i++)
		for (int j=0; j < 2*l; j++)
			{
			Edge *e = new Edge(0,i,j);
			//Edge *eref;
			//eref = &e;
			this->lattice_edges[i][j]=*e; //lattice edges holds the objects
			this->creatable.push_back(&lattice_edges[i][j]); //creatable list holds the memory location (refs)
			this->error_map[i][j]=0;

			}
	
	
	
	}

Lattice::~Lattice()
{
		delete lattice_verts;
		delete error_map;
		delete lattice_edges;

}

		
void Lattice::Add_Anyon(Anyon* a, bool u = true)
{
	int length = ArraySize;
	if (u)
	{
		anyons.push_back(a);
	}
	lattice_verts[a->x][a->y]=1;
	if (lattice_verts[mod((a->x+1),length)][a->y]==1)
	{
		lattice_edges[a->x][2*a->y].Change_Type(2); //Change lattice edge to a killable link between x,y and (x+1),y
		swapable.erase( remove( swapable.begin(), swapable.end(), &lattice_edges[a->x][2*a->y]), swapable.end() );
		killable.push_back(&lattice_edges[a->x][2*a->y]);
		lattice_edges[a->x][2*a->y].Update_Anyons(a,1);
	}
	else
	{
		lattice_edges[a->x][2*a->y].Change_Type(1); //Change lattice edge to a swapable link between x,y and (x+1),y
		creatable.erase( remove( creatable.begin(), creatable.end(), &lattice_edges[a->x][2*a->y]), creatable.end() );
		swapable.push_back(&lattice_edges[a->x][2*a->y]);
		lattice_edges[a->x][2*a->y].Update_Anyons(a,1);
	}
	
	if (lattice_verts[a->x][mod((a->y+1),length)]==1)
	{
		lattice_edges[a->x][mod(2*a->y+1,2*length)].Change_Type(2); //Change lattice edge to a killable link between x,y and x,y+1
		swapable.erase( remove( swapable.begin(), swapable.end(), &lattice_edges[a->x][mod(2*a->y+1,2*length)]), swapable.end() );
		killable.push_back(&lattice_edges[a->x][mod(2*a->y+1,2*length)]);	
		lattice_edges[a->x][mod(2*a->y+1,2*length)].Update_Anyons(a,1);
	}
	else
	{
		lattice_edges[a->x][mod(2*a->y+1,2*length)].Change_Type(1); //Change lattice edge to a swapable link between x,y and x,y+1
		creatable.erase( remove( creatable.begin(), creatable.end(), &lattice_edges[a->x][mod(2*a->y+1,2*length)]), creatable.end() );
		swapable.push_back(&lattice_edges[a->x][mod(2*a->y+1,2*length)]);
		lattice_edges[a->x][mod(2*a->y+1,2*length)].Update_Anyons(a,1);		
	}
	
	if (lattice_verts[mod((a->x-1),length)][a->y]==1)
	{
		lattice_edges[mod(a->x-1,length)][2*a->y].Change_Type(2); //Change lattice edge to a killable link between x,y and (x-1),y
		swapable.erase( remove( swapable.begin(), swapable.end(), &lattice_edges[mod(a->x-1,length)][2*a->y]), swapable.end() );
		killable.push_back(&lattice_edges[mod(a->x-1,length)][2*a->y]);
		lattice_edges[mod(a->x-1,length)][2*a->y].Update_Anyons(a,1);	
		
	}
	else
	{
		lattice_edges[mod(a->x-1,length)][2*a->y].Change_Type(1); //Change lattice edge to a swapable link between x,y and (x-1),y
		creatable.erase( remove( creatable.begin(), creatable.end(), &lattice_edges[mod(a->x-1,length)][2*a->y]), creatable.end() );
		swapable.push_back(&lattice_edges[mod(a->x-1,length)][2*a->y]);	
		lattice_edges[mod(a->x-1,length)][2*a->y].Update_Anyons(a,1);			
	}
	
	if (lattice_verts[a->x][mod((a->y-1),length)]==1)
	{
		lattice_edges[a->x][mod(2*a->y-1,2*length)].Change_Type(2); //Change lattice edge to a killable link between x,y and x,y-1
		swapable.erase( remove( swapable.begin(), swapable.end(), &lattice_edges[a->x][mod(2*a->y-1,2*length)]), swapable.end() );
		killable.push_back(&lattice_edges[a->x][mod(2*a->y-1,2*length)]);
		lattice_edges[a->x][mod(2*a->y-1,2*length)].Update_Anyons(a,1);	
	}
	else
	{
		lattice_edges[a->x][mod(2*a->y-1,2*length)].Change_Type(1); //Change lattice edge to a swapable link between x,y and x,y-1
		creatable.erase( remove( creatable.begin(), creatable.end(), &lattice_edges[a->x][mod(2*a->y-1,2*length)]), creatable.end() );
		swapable.push_back(&lattice_edges[a->x][mod(2*a->y-1,2*length)]);
		lattice_edges[a->x][mod(2*a->y-1,2*length)].Update_Anyons(a,1);	
	}
}


void Lattice::Remove_Anyon(Anyon* a, bool u = true)
{
	int length = ArraySize;
	if (u)
	{
		anyons.erase (remove( anyons.begin(), anyons.end(), a), anyons.end() );
	}
	lattice_verts[a->x][a->y]=0;
	if (lattice_verts[mod((a->x+1),length)][a->y]==1)
	{
		lattice_edges[a->x][2*a->y].Change_Type(1); //Change lattice edge to a swapable link between x,y and (x+1),y
		killable.erase( remove( killable.begin(), killable.end(), &lattice_edges[a->x][2*a->y]), killable.end() );
		swapable.push_back(&lattice_edges[a->x][2*a->y]);
		lattice_edges[a->x][2*a->y].Update_Anyons(a,0);
	}
	else
	{
		lattice_edges[a->x][2*a->y].Change_Type(0); //Change lattice edge to a killable link between x,y and (x+1),y
		swapable.erase( remove( swapable.begin(), swapable.end(), &lattice_edges[a->x][2*a->y]), swapable.end() );
		creatable.push_back(&lattice_edges[a->x][2*a->y]);
		lattice_edges[a->x][2*a->y].Update_Anyons(a,0);
	}
	
	if (lattice_verts[a->x][mod((a->y+1),length)]==1)
	{
		lattice_edges[a->x][mod(2*a->y+1,2*length)].Change_Type(1); //Change lattice edge to a swapable link between x,y and x,y+1
		killable.erase( remove( killable.begin(), killable.end(), &lattice_edges[a->x][mod(2*a->y+1,2*length)]), killable.end() );
		swapable.push_back(&lattice_edges[a->x][mod(2*a->y+1,2*length)]);
		lattice_edges[a->x][mod(2*a->y+1,2*length)].Update_Anyons(a,0);
	}
	else
	{
		lattice_edges[a->x][mod(2*a->y+1,2*length)].Change_Type(0); //Change lattice edge to a killable link between x,y and x,y+1
		swapable.erase( remove( swapable.begin(), swapable.end(), &lattice_edges[a->x][mod(2*a->y+1,2*length)]), swapable.end() );
		creatable.push_back(&lattice_edges[a->x][mod(2*a->y+1,2*length)]);
		lattice_edges[a->x][mod(2*a->y+1,2*length)].Update_Anyons(a,0);		
	}
	
	if (lattice_verts[mod((a->x-1),length)][a->y]==1)
	{
		lattice_edges[mod(a->x-1,length)][2*a->y].Change_Type(1); //Change lattice edge to a swapable link between x,y and (x-1),y
		killable.erase( remove( killable.begin(), killable.end(), &lattice_edges[mod(a->x-1,length)][2*a->y]), killable.end() );
		swapable.push_back(&lattice_edges[mod(a->x-1,length)][2*a->y]);
		lattice_edges[mod(a->x-1,length)][2*a->y].Update_Anyons(a,0);		
	}
	else
	{
		lattice_edges[mod(a->x-1,length)][2*a->y].Change_Type(0); //Change lattice edge to a killable link between x,y and (x-1),y
		swapable.erase( remove( swapable.begin(), swapable.end(), &lattice_edges[mod(a->x-1,length)][2*a->y]), swapable.end() );
		creatable.push_back(&lattice_edges[mod(a->x-1,length)][2*a->y]);
		lattice_edges[mod(a->x-1,length)][2*a->y].Update_Anyons(a,0);		
	}
	
	if (lattice_verts[a->x][mod((a->y-1),length)]==1)
	{
		lattice_edges[a->x][mod(2*a->y-1,2*length)].Change_Type(1); //Change lattice edge to a swapable link between x,y and x,y-1
		killable.erase( remove( killable.begin(), killable.end(), &lattice_edges[a->x][mod(2*a->y-1,2*length)]), killable.end() );
		swapable.push_back(&lattice_edges[a->x][mod(2*a->y-1,2*length)]);
		lattice_edges[a->x][mod(2*a->y-1,2*length)].Update_Anyons(a,0);		
	}
	else
	{
		lattice_edges[a->x][mod(2*a->y-1,2*length)].Change_Type(0); //Change lattice edge to a killable link between x,y and x,y-1
		swapable.erase( remove( swapable.begin(), swapable.end(), &lattice_edges[a->x][mod(2*a->y-1,2*length)]), swapable.end() );
		creatable.push_back(&lattice_edges[a->x][mod(2*a->y-1,2*length)]);
		lattice_edges[a->x][mod(2*a->y-1,2*length)].Update_Anyons(a,0);		
	}
}




void Lattice::Translate_Anyon(Anyon* a, int newx, int newy)
{
	//thoughts: I think I can make this function by just using Add_Anyon and Remove_Anyon cleverly.
	//All that's necessary, *I think*, is versions of those functions which don't populate/depopulate the anyons list itself
	//because I really just want to modify the edge properties
	Remove_Anyon(a,false);
	a->x = newx;
	a->y = newy;
	Add_Anyon(a,false);


}

void Lattice::Update_Sector()
{
	int xsum = 0;
	int ysum = 0;
	
	for (int i = 0; i < ArraySize; i++)
	{
		xsum += error_map[0][2*i];
		ysum += error_map[i][1];
	}
	
	SectorX = xsum%2;
	SectorY = ysum%2;
}		
		

void Lattice::Update()
{
	float Norm = creatable.size()*cRate + killable.size()*kRate + swapable.size()*sRate;
	
	float pc = creatable.size()*cRate / Norm;
	float pk = killable.size()*kRate / Norm;
//	float ps = swapable.size()*sRate / Norm;
	
	float r = rand() / (RAND_MAX + 1.);
	
	/*cout << "pc: " << pc << " pk: " << pk << " rand: " << r << endl;
	cout << "cRate: " << cRate << " kRate: " << kRate << " srate: " << sRate << endl;
	cout << "creatable.size(): " << creatable.size() << " killable.size(): " << killable.size() << " swapable.size(): " << swapable.size() << endl;
	cout << endl;*/
	
	float DeltaT = (-1./Norm) * log(r);
	
	TotalTime += DeltaT;
	
	/*if (pk > 0)
	{
		cout << "ANNIHILATION POSSIBLE!" << endl;
	}*/
	
	if (r >= 0 && r <pc)
	{
		int rce = rand() % creatable.size();	
		int xpos = creatable[rce]->x;
		int ypos = creatable[rce]->y12;
		error_map[xpos][ypos] = (error_map[xpos][ypos]+1)%2;
		
		if (ypos % 2 == 0)
		{
			Anyon* a1 = new Anyon(xpos,ypos/2);
			Anyon* a2 = new Anyon(mod(xpos+1,ArraySize),ypos/2);
			Add_Anyon(a1);
			Add_Anyon(a2);
			
		}
		else
		{
			Anyon* a1 = new Anyon(xpos,(ypos-1)/2);
			Anyon* a2 = new Anyon(xpos,mod((ypos-1)/2+1,ArraySize));
			Add_Anyon(a1);
			Add_Anyon(a2);
		}
		

	}

	if (r >= pc && r <pc+pk)
	{
		int rke = rand() % killable.size();	
		Anyon* a1 = killable[rke]->A1;
		Anyon* a2 = killable[rke]->A2;

		error_map[killable[rke]->x][killable[rke]->y12] = (error_map[killable[rke]->x][killable[rke]->y12]+1)%2;
		
		Remove_Anyon(a2);
		Remove_Anyon(a1);
		

	
	
	}
	
	if (r >= pc+pk && r <=1)
	{
		int rte = rand() % swapable.size(); //r_t_e = random translateable edge
		int xdisp = (swapable[rte]->x - swapable[rte]->A1->x);
		int ydisp = (swapable[rte]->y12 - 2*swapable[rte]->A1->y);
		
		if (xdisp == 0 && ydisp == 0)
			xdisp = 1;

		error_map[swapable[rte]->x][swapable[rte]->y12] = (error_map[swapable[rte]->x][swapable[rte]->y12] + 1)%2;
			
		Translate_Anyon(swapable[rte]->A1,mod(swapable[rte]->A1->x + xdisp,ArraySize),mod(swapable[rte]->A1->y + ydisp,ArraySize));
	

	
	}
	

	


}









void Lattice::Print_Lattice()
{
	int length = ArraySize;
	cout << endl;
	cout << "Anyon map" << endl;
	for (int j = 0; j < length; j++)
	{
		cout << endl;
		for (int i = 0; i < length; i++)
		{
			cout << lattice_verts[i][j];
		}
	}
	cout << endl;
	/*
	cout << "Edge Maps" << endl;
	
	cout << "Creatables: " << endl;
	
	for (int j = 0; j < 2*length; j++)
	{
		cout << endl;
		for (int i = 0; i < length; i++)
		{
			if (lattice_edges[i][j].edgetype == 0)
			{
				if (mod(j,2) == 0)
					cout << " _ ";
				if (mod(j,2) == 1)
					cout << "|  ";
			}
			else
			{
				if (mod(j,2) == 0)
					cout << " , ";
				if (mod(j,2) == 1)
					cout << "\\  ";
			}
		}
	}
	
	cout << endl;
	cout << "Killables: " << endl;
	
	for (int j = 0; j < 2*length; j++)
	{
		cout << endl;
		for (int i = 0; i < length; i++)
		{
			if (lattice_edges[i][j].edgetype == 2)
			{
				if (mod(j,2) == 0)
					cout << " _ ";
				if (mod(j,2) == 1)
					cout << "|  ";
			}
			else
			{
				if (mod(j,2) == 0)
					cout << " , ";
				if (mod(j,2) == 1)
					cout << "\\  ";
			}
		}
	}
	
	cout << endl;
	cout << "Swapables: " << endl;
	
	for (int j = 0; j < 2*length; j++)
	{
		cout << endl;
		for (int i = 0; i < length; i++)
		{
			if (lattice_edges[i][j].edgetype == 1)
			{
				if (mod(j,2) == 0)
					cout << " _ ";
				if (mod(j,2) == 1)
					cout << "|  ";
			}
			else
			{
				if (mod(j,2) == 0)
					cout << " , ";
				if (mod(j,2) == 1)
					cout << "\\  ";
			}
		}
	}
	*/
	cout << endl;
	cout << "List data:" << endl;
	cout << "csize:" << creatable.size() << " ssize:" << swapable.size() << " ksize:" << killable.size() << endl;
	cout << "anyons: " << anyons.size() << endl;
	cout << "Sector: " << SectorX << SectorY << endl;
	cout << endl;
}


	
		



