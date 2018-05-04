/*  Copyright (C) 2011 Jens Pfau <jpfau@unimelb.edu.au>
 *  Copyright (C) 2014 Alex Stivala <stivalaa@unimelb.edu.au>
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <string>
#include <limits>
#include <cassert>
#include <unistd.h> // for getpid()


#include "model.hpp"

#define SVN_VERSION_MODEL_CPP "lattice dyadic influence von Neumann model with noise rev 16 $Id: model.cpp 783 2016-07-03 08:18:47Z stivalaa $"

// this version has the social network (graph) removed for
// efficiency on simpler model that doesnot use it

// modified by ADS to seed srand with time+pid and to also 
// ensure unique tmp files so safe for parallel execution with MPI Python
// NB this involved extra command line pameter to model to not compatible
// with unmodified versions
// Also added theta parameter as threshold for bounded confidence model



static char *tempfileprefix = NULL;


// Read the time_list---the iterations at which the stats of the simulation are
// to be printed out---from a temporary file.
int read_time_list(unsigned long long int** time_list, int* n) {
        std::ifstream inFile ((toString(tempfileprefix) + ".T").c_str());

	if (inFile) {
		std::string line;

		// Get the first line and read out the number of time steps in the list
		if (!getline(inFile, line))
			return -1;

		*n = convert<unsigned long long int>(line);
		unsigned long long int* tmplist = new unsigned long long int[*n];

		// Get the list itself
		if (!getline(inFile, line))
			return -1;

		// Read every item of the list
		std::istringstream linestream(line);
		std::string item;
		int i = 0;
		while (getline (linestream, item, ',') and i < *n) {
			tmplist[i] = convert<unsigned long long int>(item);
			i++;
		}

		if (i != *n)
			return -1;

		*time_list = tmplist;
	} else {
		*n = 0;
	}

	return 0;
}



// Calculate Manhattan distance on the lattice between locatino (x1,y1) 
// and (x2,y2)
inline int manhattan_distance(int x1, int y1, int x2, int y2) {
  return abs(x2 - x1) + abs(y2 - y1);
}



// Calculate the cultural similarity between agents a and b as the proportion of
// features they have in common.
inline double similarity(Grid<int>& C, const int F, const int a, const int b) {
	int same = 0;
	for (int i = 0; i < F; i++)
			same += (C(a,i) == C(b,i));
	return (double)same/F;
}



// test for equilbrium, return true if no more change is possible
bool stop2(Grid<int>& C, Grid<int>& L, int n, int F, double theta,
           int vonNeumann_radius) {
    // at equilibrium, all agents in a neighbourhood
    // must have either identical culture, or
    // or completely distinct culture (no traits in common, similarity = 0), or
    // cultures with similarity < theta.
    assert(theta >= 0.0 && theta <= 1.0);
    for (int i = 0; i < n; i++)  {
        for (int j = i+1; j < n; j++) { // symmetric, only need i,j where j > i
            if (manhattan_distance(L(i,0), L(i,1), L(j, 0), L(j,1)) > vonNeumann_radius) {
              // TODO this could be made more efficient by looping over only
              // agents in neighbourhood instead of this way of looping over
              // all and just doing next if not in neighbourhood.
              // But doesn't much matter in this versino since we only
              // use this when noise=0, in all other cases have to go
              // to iteration limit tmax anyway.
              continue;
            }
            double sim = similarity(C, F, i, j);
            assert(sim >= 0.0 && sim <= 1.0);
            // NB using >= and <= not == to for floating point comparison
            // with 0 or 1 since floating point == is dangerous, but
            // values cannot be < 0 or > 1, as just asserted so equivalent
            // to equality
            if ( !((sim >= 1.0 ) ||
                   (sim <= 0.0 || sim < theta )) ) {
                return false;
            }
        }
    }
    return true;
}



//
// the model main loop
// 
unsigned long long model(Grid<int>& L, Grid<int>& O, Grid<int>& C,
                         Grid<int>&Strategy, // not actually used
                         unsigned long long int tmax, int n, int m, int F,
		int q, double r, double s,
		bool toroidal, bool network, double tolerance, bool directed_migration,
		double phy_mob_a, double phy_mob_b, double soc_mob_a, double soc_mob_b,
		double r_2, double s_2, double phy_mob_a_2, double phy_mob_b_2, double soc_mob_a_2, double soc_mob_b_2,
		double k, unsigned long long int timesteps[], long int time_list_length, std::ofstream& log, 
                         double theta, int vonNeumann_radius, double noise_level
  ) {
       srand(time(NULL)+(unsigned int)getpid()); // so processes started at same time have different seeds (time is only second resolution)



	std::cout << "tmax: " << tmax << std::endl;
	std::cout << "n: " << n << std::endl;
	std::cout << "m: " << m << std::endl;
	std::cout << "F: " << F << std::endl;
	std::cout << "q: " << q << std::endl;
  std::cout << "theta: " <<  theta << std::endl;
  std::cout << "radius: " <<  vonNeumann_radius << std::endl;
  std::cout << "noise: " << noise_level << std::endl;


    int a, b, idx;
    int nextstep = 0;

  // at the moment always have r = r' as in Flache & Macy (2011)
  double interaction_noise_level = noise_level; // r in Flache & Macy (2011)



	// run model
    for (unsigned long long t = 0; t < tmax; t++) {
    	// If this iteration is in the time list, write out the current state
    	// of the simulation.
    	if (nextstep < time_list_length and timesteps[nextstep] == t) {
        L.write((toString(tempfileprefix) + "-" + toString(timesteps[nextstep]) + ".L").c_str(), ',');
        C.write((toString(tempfileprefix) + "-" + toString(timesteps[nextstep]) + ".C").c_str(), ',');
        Strategy.write((toString(tempfileprefix) + "-" + toString(timesteps[nextstep]) + ".Strategy").c_str(), ','); // not used but stops python script crashing
    		nextstep++;
    	}
      

    	if (t == 50000 || t == 100000 || t == 500000 || t == 1000000 || (t > 0 && t % 10000000 == 0)) {
    		std::cout << "Reaching " << t << " iterations." << std::endl;
        // Only check absorbing state if there is no noise and not doing time series
			  if (noise_level == 0.0 && time_list_length == 0 && stop2(C, L, n, F, theta, vonNeumann_radius)) {
  		  		std::cout << "Stopping after " << t << " iterations." << std::endl;
            return t;
  			}
      }


    	// Draw one agent randomly.
    	a = rand() % n;

      // Make lsit of neigbhours in von Neumann
      // neighbourhood of the focal agent.
      int ax = L(a,0);
      int ay = L(a,1);
      std::vector<int>neighbours;
      for (int xi = -1*vonNeumann_radius; xi <= vonNeumann_radius; xi++) {
          for (int yi = -1*vonNeumann_radius; yi <= vonNeumann_radius; yi++) {
              if ((xi != 0 || yi != 0) &&  // don't include focal agent itself
                  abs(xi) + abs(yi) <= vonNeumann_radius) {
                  int bx = ax + xi;
                  int by = ay + yi;
                  // handle edge of lattice, not toroidal
                  if (bx >= 0 && bx < m && by >= 0 && by < m) {
                      neighbours.push_back(O(bx, by));
                  }
              }
          }
      }
     
      // choose a random agent 
      // If not identical to focal agent then
      // randomly decide on one feature that a and b do not have in common yet.
      b = neighbours[rand() % neighbours.size()];
      double sim = similarity(C, F, a, b);
      // if sim is at least theta (bounded confidence)
      // a and b interact with probability sim
      if (sim >= theta && rand()/(float)RAND_MAX < sim) {
        if (sim < 1.0) {
          do {
            idx = rand() % F;
          } while (C(a,idx) == C(b, idx));
          // Let a copy this feature from b.
          C(a,idx) = C(b,idx);

        }
      }
      // as in KETS step 6 in Flache & Macy (2011) interaction noise
      // is introduced by changing a random trait to a random value
      // with probability interaction_noise_level (r)
      if (rand()/(float)RAND_MAX < interaction_noise_level) {
        idx = rand() % F;
        C(a, idx) = rand() % q;
      }
    }
    std::cout << "Stopping after tmax = " << tmax << " iterations." << std::endl;
    return tmax;
}




int main(int argc, char* argv[]) {
	std::ofstream log("log.txt");

	// If the binary file is called with the argument -v, only the svn version
	// this binary was compiled from is printed.
	if (argc == 2 and argv[1][0] == '-' and argv[1][1] == 'v') {
		std::cout << "model.hpp: " << SVN_VERSION_MODEL_HPP << ", model.cpp: " << SVN_VERSION_MODEL_CPP << std::endl;
		return 0;
	}


	// Otherwise set default model arguments.
	int n = 100, m = 10, F = 5, q = 15;

	unsigned long long int tmax = 100000;
	double r = 1;
	double s = 1;

	bool toroidal = false;
	bool network = false;

    double tolerance = -1;
    bool directed_migration = false;

    double phy_mob_a = 1;
    double phy_mob_b = 10;
    double soc_mob_a = 1;
    double soc_mob_b = 10;

	double r_2 = r;
	double s_2 = s;

    double phy_mob_a_2 = phy_mob_a;
    double phy_mob_b_2 = phy_mob_b;
    double soc_mob_a_2 = soc_mob_a;
    double soc_mob_b_2 = soc_mob_b;

    double k = 0.01;

   double theta = 0.0;
    // size of the von Neumann neighbourhood for actors to interact, i.e.
    // the maximum Manhattan distance between interacting actors
    int vonNeumann_radius = 1;

    double noise_level = 0.0; // rate of noise

    int num_joint_activities = 10; // number of different public goods games
    int pool_multiplier = 1; // multiplier of pool value in public goods games

    

    // If there are arguments, assume they hold model arguments in the following
    // order.
	if (argc > 1) {
		int index = 1;
		tmax = atoll(argv[index++]);
		n = atoi(argv[index++]);
		m = atoi(argv[index++]);
		F = atoi(argv[index++]);
		q = atoi(argv[index++]);
		r = atof(argv[index++]);                     // not used
		s = atof(argv[index++]);                     // not used
		toroidal = atoi(argv[index++]);              // not used
		network = atoi(argv[index++]);               // not used
		tolerance = atof(argv[index++]);             // not used
		directed_migration = atoi(argv[index++]);    // not used
		phy_mob_a = atof(argv[index++]);             // not used
		phy_mob_b = atof(argv[index++]);             // not used
		soc_mob_a = atof(argv[index++]);             // not used
		soc_mob_b = atof(argv[index++]);             // not used

		r_2 = atof(argv[index++]);                   // not used
		s_2 = atof(argv[index++]);                   // not used
		phy_mob_a_2 = atof(argv[index++]);           // not used
		phy_mob_b_2 = atof(argv[index++]);           // not used
		soc_mob_a_2 = atof(argv[index++]);           // not used
		soc_mob_b_2 = atof(argv[index++]);           // not used
		k = atof(argv[index++]);                     // not used
		tempfileprefix = argv[index++];
    theta = atof(argv[index++]);
    vonNeumann_radius = atoi(argv[index++]);   // not used: von Neumann radius
    noise_level = atof(argv[index++]);         // noise level
    num_joint_activities = atoi(argv[index++]);    // not useed:  number of games
    pool_multiplier = atoi(argv[index++]); // not used:  pool multiplier
	}


    if (toroidal) {
    	std::cout << "alarm, toroidal not supported at the moment" << std::endl;
    	return -1;
    }

   if (n != m*m) {
    std::cerr << "for von Neumann neighbourrhood model must have n =m*m so all latice points used, but m = " << m << " and n = " << n <<std::endl;
    return -1;
   }
	// Try to read the list of iterations that determine when statistics are to
	// be created from a temporary file.
	unsigned long long int* time_list = NULL;
	int time_list_length = 0;
	int res = read_time_list(&time_list, &time_list_length);
	if (res == -1) {
		std::cout << "The time list file could not be read or there was a problem with its format." << std::endl;
		return -1;
	}

	Grid<int> L(n,2,-1);
	Grid<int> C(n,F,0);
	Grid<int> O(m,m,-1); // in this version O(x,y) contains agent id at that location or -1 if unoccupied
        Grid<int>Strategy(n,1,STRATEGY_COOPERATE); // not used but stops python scripts crashing



	if (argc == 1) {
		std::cerr << "must provide initialization parameters" << std::endl;
	} else {
		// load data from file
          //G.read((toString(tempfileprefix) + ".adj").c_str(), ' ');
 	        L.read((toString(tempfileprefix) + ".L").c_str(), ',');
	        C.read((toString(tempfileprefix) + ".C").c_str(), ',');
                Strategy.read((toString(tempfileprefix) + ".Strategy").c_str(), ',');; // not used but stops python scripts crashing
	}

	for (int i = 0; i < n; i++)
		O(L(i,0), L(i,1)) = i; // O(x,y) is agent at lattice location (x,y)


	// Call the model
        unsigned long long tend = model(L, O, C, Strategy, tmax, n, m, F, q, r, s, toroidal, network, tolerance, directed_migration,
    		phy_mob_a, phy_mob_b, soc_mob_a, soc_mob_b,
    		r_2, s_2, phy_mob_a_2, phy_mob_b_2, soc_mob_a_2, soc_mob_b_2, k,
    		time_list, time_list_length, log, theta, vonNeumann_radius,
                noise_level);  

        std::cout << "Last iteration: " << tend << std::endl;
    
        // Write out the state of the simulation
        L.write((toString(tempfileprefix) + ".L").c_str(), ',');
        C.write((toString(tempfileprefix) + ".C").c_str(), ',');
        Strategy.write((toString(tempfileprefix) + ".Strategy").c_str(), ','); // not used but stops python scripts crashing

	std::ofstream outFile((toString(tempfileprefix) + ".Tend").c_str());
	outFile << tend;


	delete[] time_list;
	time_list = NULL;

	std::cout << "Fin" << std::endl;
	return 0;
}
