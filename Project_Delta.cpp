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
void print(vector<ship> fleet_dat);


//===========================================================================					main
int main()
{
	// Seed Random Number Generator
	srand(time(NULL));

	// Welcome Screen
	welcome();

	// Set Timestep Resolution
	double dt = .2;

	// Set Maximum Simulation Time (Seconds)
	double max_time = 150;

	// Number of Generations
	int gen_num = 10;

	// Initialize Ships
	int ship_num = 5;
	fleet starfleet;

	for (int i = 0; i < 2*ship_num; i++)
	{
		ship drydock;
		drydock.set_inital_conditions();
		starfleet.fleet_dat.push_back(drydock);
	}

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

		starfleet.downselect();
		cout << size(starfleet.fleet_dat) << '\t';

		starfleet.repopulate();
		cout << size(starfleet.fleet_dat) << endl;
	}

	// Print to File
	print(starfleet.fleet_dat);

    return 0;
}

/// Main Functions
//===========================================================================					welcome
void welcome()
{
	cout << "Ship Simulator Running\n\n" << endl;
}

//===========================================================================					print
void print(vector<ship> fleet_dat)
{
	// Clear Console Screen
	// system("CLS");

	// Now Writing Data To File
	cout << "Now Writing Data To File" << endl;

	// Create output file
	ofstream output_file;
	output_file.open("Captains_Log.txt");

	int end = fleet_dat.size();
	for (int i = 0; i < end; i++)
	{
		cout << i << '\t' << fleet_dat.at(i).fitness << endl;
		//output_file << i << '\t' << fleet_dat.at(i).fitness << endl;
	}

	int j = 0;
	end = fleet_dat.at(j).x.size();
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

	//User console update
	cout << "Data Has Been Written To File" << endl;
}

/// Fleet Functions
//===========================================================================					repopulate
void fleet::repopulate()
{
	size_t L = size(fleet_dat);
	double Len = L;

	for (int i = 0; i < Len; i++)
	{
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
	x.push_back(50);
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
		double we = ((rand() % 100)) - ((rand() % 100));
		weights.push_back(we);
	}
	// Set Compass
	calc_compass_heading();
}

//===========================================================================					mutate
void ship::mutate()
{
	fitness = 0;

	for (int i = 0; i < weight_num; i++)
	{
		weights.at(i) = weights.at(i) + (((rand() % 100)) - ((rand() % 100)));
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
	cout << compass << endl;
}

//===========================================================================					simulate
void ship::simulate(double dt, double max_time)
{
	// Set Up Neural Network
	neural_network NN;
	NN.setup(1, 5, 1);
	NN.set_in_min_max(-180, 180);
	NN.set_out_min_max(-26, 26);
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
		u = NN.get_output(0) *.01;
		omega = (omega + (u - omega * dt / T));
		theta.push_back(theta.at(i) + omega * dt);
		y.push_back(y.at(i) + speed * sin(theta.at(i + 1)) * dt);
		x.push_back(x.at(i) + speed * cos(theta.at(i + 1)) * dt);
		i++;
		stime = stime + dt;
		
		// Fitness of the Ships Captain
		calc_compass_heading();
		double compass_r = (compass * 3.1415) / 180;
		fitness = fitness + cos(compass_r) - 1;

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
						fitness = fitness + 200;
						stime = max_time;
						cout << "Goal Reached" << endl;
						break;
					}
				}
			}
		}

		// Check for Out of Bounds Ship
		if (x.at(i) > 1000 || x.at(i) < 0 || y.at(i) > 1000 || y.at(i) < 0)
		{
			fitness = fitness - 5000;
			stime = max_time;
			break;
		}
	}

	stime = 0;
}






