//===========================================================================
//Delta.cpp
//
//
//Input:
//	N/A
//			
//Output:
//  Captains_Log.txt
//
//===========================================================================

#include "stdafx.h"
#include <conio.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <cassert>
#include <random>
#include <vector>
#include <iomanip>
#include <deque>
#include <algorithm>
#include <numeric>
#include <functional>
#include <Windows.h>
#include <array>
#include "LY_NN.h"

using namespace std;

// Classes
class ship
{
public:
	double fitness;
	double compass;
	double speed;
	double omega;
	double weight_num;
	double goal[3];
	double stime;
	vector<double> x;
	vector<double> y;
	vector<double> theta;
	vector<double> weights;

	void set_inital_conditions();
	void mutate();
	void calc_compass_heading();
	void simulate(double dt, double max_time);
};
class fleet
{
public:
	vector<ship> fleet_dat;
	void repopulate();
	void downselect();
};

// Function Prototypes
void welcome();
void print(vector<ship> fleet_dat, vector<ship> every_five);


//===========================================================================					main
int main()
{
	// Seed Random Number Generator
	srand(time(NULL));

	// Welcome Screen
	welcome();

	// Set Timestep Resolution
	double dt = .4;

	// Set Maximum Simulation Time (Seconds)
	double max_time = 100;

	// Number of Generations
	int gen_num = 20;

	// Initialize Ships
	int ship_num = 20;
	fleet starfleet;

	for (int i = 0; i < 2*ship_num; i++)
	{
		ship drydock;
		drydock.set_inital_conditions();
		starfleet.fleet_dat.push_back(drydock);
	}
	vector<ship> every_five;

	// Loop through Generations
	for (int i = 0; i < gen_num; i++)
	{
		cout << "Gen Number " << i + 1 << endl;
		// Simulate Ships and Set Fitness
		int j = 0;
		for (j = 0; j < 2 * ship_num; j++)
		{
			//cout << j << endl;
			starfleet.fleet_dat.at(j).simulate(dt, max_time);
		}
		if (i % 5 == 0)
		{
			every_five.push_back(starfleet.fleet_dat.at(j));
		}
		if (i + 1 < gen_num)
		{
		starfleet.downselect();
		cout << size(starfleet.fleet_dat) << '\t';

		starfleet.repopulate();
		cout << size(starfleet.fleet_dat) << endl;
		}
		cout << 'test' << endl;
		if (i % 5 == 0)
		{
			int end = starfleet.fleet_dat.size();
			int best = 0;
			for (int i = 0; i < end; i++)
			{
				if (i > 0)
				{
					if (starfleet.fleet_dat.at(i).fitness >= starfleet.fleet_dat.at(best).fitness)
					{
						best = i;
					}
				}
			}
			every_five.push_back(starfleet.fleet_dat.at(best));
		}
	}

	// Print to File
	print(starfleet.fleet_dat, every_five);

    return 0;
}

/// Main Functions
//===========================================================================					welcome
void welcome()
{
	cout << "Ship Simulator Running\n\n" << endl;
}

//===========================================================================					print
void print(vector<ship> fleet_dat, vector<ship> every_five)
{
	// Clear Console Screen
	// system("CLS");

	// Now Writing Data To File
	cout << "Now Writing Data To File" << endl;

	// Create output file
	ofstream output_file;
	output_file.open("Captains_Log.txt");

	int end = fleet_dat.size();
	int best = 0;
	for (int i = 0; i < end; i++)
	{
		cout << i << '\t' << fleet_dat.at(i).fitness << endl;
		//output_file << i << '\t' << fleet_dat.at(i).fitness << endl;
		if (i > 0)
		{
			if (fleet_dat.at(i).fitness >= fleet_dat.at(best).fitness)
			{
				best = i;
			}
		}
	}

	int j = best;
	end = fleet_dat.at(j).x.size();
	output_file << j << '\t' << 0000 << '\t' << endl;
	for (int i = 0; i < end; i++)
	{
		output_file << fleet_dat.at(j).x.at(i) << '\t' << fleet_dat.at(j).y.at(i) << endl;
	}

	//Close output file
	output_file.close();

	// Wait for
	//Sleep(750);
	//system("PAUSE");

	// Clear Console Screen
	//system("CLS");

	// Create output file
	ofstream out_file;
	out_file.open("Captains_Log_Full.txt");

	end = every_five.at(j).x.size();
	int walk = every_five.size();
	for (int i = 0; i < end; i++)
	{
		for (int k = 0; k < walk; k++)
		{
			out_file << fleet_dat.at(k).x.at(i) << '\t' << fleet_dat.at(k).y.at(i) << '\t\t';
		}
		out_file << endl;
	}

	//Close output file
	out_file.close();

	//User console update
	cout << "Data Has Been Written To File" << endl;
}

/// Fleet Functions
//===========================================================================					repopulate
void fleet::repopulate()
{
	size_t L = size(fleet_dat);
	double Len = L;
	cout << "repo" << endl;
	for (int i = 0; i < Len; i++)
	{
		fleet_dat.at(i).fitness = 0;
		ship spare = fleet_dat.at(i);
		spare.mutate();
		fleet_dat.push_back(spare);
	}
	assert(size(fleet_dat) == 2*L);
}

//===========================================================================					downselect
void fleet::downselect()
{
	size_t L = size(fleet_dat);
	double Len = L;
	
	for (int i = 0; i < .5 * Len; i++)
	{
		int opponent_1 = rand() % size(fleet_dat);
		int opponent_2 = rand() % size(fleet_dat);

		while (opponent_1 == opponent_2)
		{
			if (opponent_2 != size(fleet_dat) - 1)
			{
				opponent_2 = opponent_2 + 1;
			}
			else
			{
				opponent_2 = opponent_2 - 1;
			}
		}

		if (fleet_dat.at(opponent_1).fitness >= fleet_dat.at(opponent_2).fitness)
		{
			fleet_dat.erase(fleet_dat.begin() + opponent_2);

		}
		else
		{
			fleet_dat.erase(fleet_dat.begin() + opponent_1);
		}
	}

	assert(size(fleet_dat) == .5*L);
	
}

/// Ship Functions
//===========================================================================					set_inital_conditions
void ship::set_inital_conditions()
{
	// Set Initial stime
	stime = 0;

	// Initial Fitness
	fitness = 0;

	// Initial Position of ship
	x.push_back(800);
	y.push_back(500);


	// Position of Goal
	goal[0] = 900;
	goal[1] = 500;
	goal[2] = 50;

	// Initial Speed
	speed = 3;

	// Initial Direction & Angular Velocity
	theta.push_back(0);
	omega = 0;

	// Set Weights
	weight_num = 16;

	for (int i = 0; i < weight_num; i++)
	{
		double we = ((rand() % 1000)) - ((rand() % 1000));
		weights.push_back(we);
	}
	// Set Compass
	calc_compass_heading();
}

//===========================================================================					mutate
void ship::mutate()
{
	fitness = 0;
	x.clear();
	y.clear();
	theta.clear();

	// Initial Position of ship
	x.push_back(800);
	y.push_back(500);
	theta.push_back(0);

	for (int i = 0; i < weight_num; i++)
	{
		weights.at(i) = weights.at(i) + (((rand() % 1000)) - ((rand() % 1000)));
	}
}

//===========================================================================					calc_compass_heading
void ship::calc_compass_heading()
{
	//cout << (goal[1] - y.back()) << '\t' << (goal[0] - x.back()) << '\t';
	compass = atan((goal[1] - y.back()) / (goal[0] - x.back()));
	//cout << theta.back() << endl;
	compass = (compass * (double)180) / 3.1415;
	int it = 0;
	///*
	double com = theta.back();
	while (it != 1)
	{
		if (com >= 180)
		{
			com = com - 180;
		}
		if (com <= -180)
		{
			com = com + 180;
		}
		if (com >= -180 && com <= 180)
		{
			it = 1;
		}
	}
	//*/
	compass = com - compass;
	//cout << compass << endl;
}

//===========================================================================					simulate
void ship::simulate(double dt, double max_time)
{
	// Set Up Neural Network
	neural_network NN;
	NN.setup(1, 5, 1);
	NN.set_in_min_max(-180, 180);
	NN.set_out_min_max(-15, 15);
	NN.set_weights(weights, true);

	vector<double> compass_temp;
	compass_temp.push_back(0);
	int i = 0;
	double T = 5;
	double u = .0000;
	while (stime < max_time)
	{

		// Euler's to the Next Timestep
		compass_temp.at(0) = compass;
		//cout << compass << endl;
		NN.set_vector_input(compass_temp);
		NN.execute();
		u = NN.get_output(0) *.001 ;
		omega = (omega + (u - omega * dt / T));
		theta.push_back(theta.at(i) + omega * dt);
		y.push_back(y.at(i) + speed * sin(theta.at(i + 1)) * dt);
		x.push_back(x.at(i) + speed * cos(theta.at(i + 1)) * dt);
		i++;
		stime = stime + dt;
		
		// Fitness of the Ships Captain
		calc_compass_heading();
		double compass_r = (compass * 3.1415) / 180;
		fitness = fitness + (cos(compass_r) - 1) * .1;

		// Check for Goal Passage
		double distance_in_one_step = speed * dt * 2;

		if (x.at(i) >= goal[0] - distance_in_one_step)
		{
			if (x.at(i) < goal[0] + distance_in_one_step)
			{
				if (y.at(i) >= goal[1] - goal[2])
				{
					if (y.at(i) < goal[1] + goal[2])
					{
						fitness = fitness + 20000;
						stime = max_time + 1;
						cout << "Goal Reached " << fitness << endl;
					}
				}
			}
		}

		// Check for Out of Bounds Ship
		if (x.at(i) > 1000 || x.at(i) < 0 || y.at(i) > 1000 || y.at(i) < 0)
		{
			fitness = fitness - 5000;
			stime = max_time + 1;
			break;
		}
	}

	stime = 0;
}






