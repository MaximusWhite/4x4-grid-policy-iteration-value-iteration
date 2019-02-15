// ass2.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include <iostream>
#include <algorithm> 
#include <vector>
#include <cmath>

using namespace std;
//consts

enum Direction {UP, DOWN, LEFT, RIGHT, IDLE, OFF};
const double p1 = 0.5;
const double p2 = 0.3;
const double p3 = (1 - p1 - p2) / 2;
const double off_the_grid_ignore = p1 + p2;
const double off_the_grid_adjacent = (1 - p1 - p2) / 2;
const double gamma = 0.95;
const double theta = 0.001;
const vector<double> rewards = { 1.0, 2.0, 5.0, 3.1, 0.0 };

vector<double> V(16, 0.0); // initial V(s) = 0, for all s
vector<vector<double>> initial_policies(16, vector<double>(4, 0.25)); // policy is equiprobable for each action in any given state 

int get_direction(int state, int state_prime) {
	if (state_prime == state - 4) return UP;
	if (state_prime == state + 4) return DOWN;
	if (state_prime == state + 1) return RIGHT;
	if (state_prime == state - 1) return LEFT;
	if (state_prime == state) return IDLE;
	return -1;
}

bool opposite_directions(int dir1, int dir2) {
	return (dir1 == UP && dir2 == DOWN) ||
		(dir1 == DOWN && dir2 == UP) ||
		(dir1 == LEFT && dir2 == RIGHT) ||
		(dir1 == RIGHT && dir2 == LEFT);
}

double action_probability(int state, int state_prime, int action) {
	int move_direction = get_direction(state, state_prime);
	if (move_direction == IDLE) return p2 * off_the_grid_ignore;
	if (action == move_direction) {
		return p1;
	} else {
		if (move_direction != -1 && !opposite_directions(move_direction, action)) {
			return p2 * p3; 
		}
	}
	return 0.0;
}

void policy_evaluation() {
	double delta;
	int iterations = 0;
	do {
		delta = 0;
		for (int s = 0; s < 16; s++) {  // for all s 
			double v = V[s];
			double sum1 = 0.0;
			for (int a = 0; a < 4; a++) { // for all actions
				double policy_a_given_s = initial_policies[s][a];

				double sum2 = 0.0;
				for (int s_prime = 0; s_prime < 16; s_prime++) { // for all next states

					double prob = action_probability(s, s_prime, a);
					int direction = get_direction(s, s_prime);
					double reward = direction == -1 ? 0.0 : rewards[direction];
					sum2 += (prob * (reward + gamma * V[s_prime]));
				}
				sum1 += (policy_a_given_s * sum2);
			}

			V[s] = sum1;
			delta = max(delta, abs(v - V[s]));
		}

		iterations++;
		system("CLS");
		cout << "Iterations: " << iterations << endl;
	} while (delta > theta);
}

bool policy_improvement() {
	bool stable = true;
	
	for (int s = 0; s < 16; s++) {
		int max = 0;
		double max_val = 0.0;
		for (int a = 0; a < 4; a++) {
			for (int s_prime = 0; s_prime < 16; s_prime++) {
				double prob = action_probability(s, s_prime, a);
				int direction = get_direction(s, s_prime);
				double reward = direction == -1 ? 0.0 : rewards[direction];


			}
		}

	}
}

void policy_iteration() {
	bool policy_stable = false;
	do {
		policy_evaluation();
		policy_stable = policy_improvement();
	} while (!policy_stable);
}

int main(int argc, char *argv[])
{
	policy_iteration();
	cout << "Values: " << endl;
	for (int i = 0; i < 16; i++) {
		cout << "State: " << i + 1 << "; Value: " << V[i] << endl;
	}

}

