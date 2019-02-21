#include <iostream>
#include <algorithm> 
#include <vector>
#include <cmath>
#include <ctime>

using namespace std;
//consts

enum Direction {UP, DOWN, LEFT, RIGHT, IDLE, DIAG_UL, DIAG_UR, DIAG_DL, DIAG_DR, OFF};
double p1 = 0.5;
double p2 = 0.4;
const double p3 = (1 - p1 - p2) / 2;
const double off_the_grid_ignore = p1 + p2;
const double off_the_grid_adjacent = (1 - p1 - p2) / 2;
const double gamma_const = 0.95;
const double theta = 0.001;
vector<double> rewards = { -5.0, -1.0, -1.0, -1.0};

vector<double> V(16, 0.0); // initial V(s) = 0, for all s
vector<vector<vector<double>>> transition_probabilities(16, vector<vector<double>>(4, vector<double>(16, 0.0))); // s, a, s'
vector<vector<double>> policies(16, vector<double>(4, 0.25)); // policy is equiprobable for each action in any given state 

int get_direction(int state, int state_prime);
void initialize_probabilities();
int is_on_edge(int state);
double action_probability(int state, int state_prime, int action);
void policy_evaluation();
bool policy_improvement();
void policy_iteration();
int max_index(vector<double> arr);

int max_index(vector<double> arr) {
	int max = 0;
	double max_val = arr[max];
	for (int i = 1; i < arr.size(); i++) {
		if (arr[i] > max_val) {
			max = i; 
			max_val = arr[i];
		}
	}
	return max;
}

int get_direction(int state, int state_prime) {
	if (state_prime == state - 4) return UP;
	if (state_prime == state + 4) return DOWN;
	if (state_prime == state + 1) return RIGHT;
	if (state_prime == state - 1) return LEFT;
	if (state_prime == state - 5) return DIAG_UL;
	if (state_prime == state - 3) return DIAG_UR;
	if (state_prime == state + 3) return DIAG_DL;
	if (state_prime == state + 5) return DIAG_DR;
	if (state_prime == state) return IDLE;
	return -1;
}

bool opposite_directions(int dir1, int dir2) {
	return (dir1 == UP && dir2 == DOWN) ||
		(dir1 == DOWN && dir2 == UP) ||
		(dir1 == LEFT && dir2 == RIGHT) ||
		(dir1 == RIGHT && dir2 == LEFT);
}

void initialize_probabilities() {
	for (int s = 0; s < 16; s++) {
		for (int a = 0; a < 4; a++) {
			for (int s_prime = 0; s_prime < 16; s_prime++) {
				transition_probabilities[s][a][s_prime] = action_probability(s, s_prime, a);
			}
		}
	}
}

int is_on_edge(int state) {
	if (state == 0 || state == 3) return UP;
	if (state == 13 || state == 14) return DOWN;
	if (state == 4 || state == 8) return LEFT;
	if (state == 7 || state == 11) return RIGHT;
	if (state == 0) return DIAG_UL;
	if (state == 3) return DIAG_UR;
	if (state == 12) return DIAG_DL;
	if (state == 15) return DIAG_DR;
	return -1;
}

double action_probability(int state, int state_prime, int action) {
	
	double prob = 0.0;

	if (state == 0 || state == 15) {
		return prob;
	}

	int move_direction = get_direction(state, state_prime); // position of s' relative to s
	int edge_position = is_on_edge(state);

	switch (action) {
		/* ACTION = UP ********************************************************** */
		case UP:
			switch (move_direction) {   // hypothetical move from state s to s'
				case UP:			// if s' is up from s
					if (edge_position == LEFT || edge_position == RIGHT) {
						prob = p1 + p3;
					} else {
						prob = p1;
					}
					break;
				case IDLE: { 	// idle meaning going from state s to state s again
					switch (edge_position) // depending on which edge the state is (if any), the idle probability will be different
					{
						case UP: prob = p1 + p2;	// if at the top edge, that means that agent is trying to go off the grid
								break;
						case DIAG_UL: prob = p1 + p2;
								break;
						case DIAG_UR: prob = p1 + p2;
								break;
						default: prob = p2;
								break;
					}
				}
				case LEFT: if (edge_position == UP) prob = p3;
						break;
				case RIGHT: if (edge_position == UP) prob = p3;
						break;
				case DIAG_UL: prob = p3; 
						break;
				case DIAG_UR: prob = p3; 
						break;
			};
			break;
		/* ACTION = DOWN ********************************************************** */
		case DOWN:
			switch (move_direction) {   // hypothetical move from state s to s'
				case DOWN:			// if s' is up from s
					if (edge_position == LEFT || edge_position == RIGHT) {
						prob = p1 + p3;
					} else {
						prob = p1;
					}
					break;
				case IDLE: { 	// idle meaning going from state s to state s again
					switch (edge_position) // depending on which edge the state is (if any), the idle probability will be different
					{
						case DOWN: prob = p1 + p2;	// if at the bottom edge, that means that agent was trying to go off the grid
								break;
						case DIAG_DL: prob = p1 + p2;
								break;
						case DIAG_DR: prob = p1 + p2;
								break;
						default: prob = p2;
								break;
					}
				}
				case LEFT: if (edge_position == DOWN) prob = p3;
						break;
				case RIGHT: if (edge_position == DOWN) prob = p3;
						break;
				case DIAG_UL: prob = p3; 
						break;
				case DIAG_UR: prob = p3; 
						break;
			};
			break;

		/* ACTION = LEFT ********************************************************** */
		case LEFT:
			switch (move_direction) {   // hypothetical move from state s to s'
				case LEFT:			// if s' is up from s
					if (edge_position == UP || edge_position == DOWN) {
						prob = p1 + p3;
					} else {
						prob = p1;
					}
					break;
				case IDLE: { 	// idle meaning going from state s to state s again
					switch (edge_position) // depending on which edge the state is (if any), the idle probability will be different
					{
						case LEFT: prob = p1 + p2;	// if at the left edge, that means that agent is trying to go off the grid
								break;
						case DIAG_UL: prob = p1 + p2;
								break;
						case DIAG_DL: prob = p1 + p2;
								break;
						default: prob = p2;
								break;
					}
				}
				case UP: if (edge_position == LEFT) prob = p3;
						break;
				case DOWN: if (edge_position == LEFT) prob = p3;
						break;
				case DIAG_UL: prob = p3; 
						break;
				case DIAG_DL: prob = p3; 
						break;
			};
			break;

		/* ACTION = RIGHT ********************************************************** */
		case RIGHT:
			switch (move_direction) {   // hypothetical move from state s to s'
				case RIGHT:			// if s' is up from s
					if (edge_position == UP || edge_position == DOWN) {
						prob = p1 + p3;
					} else {
						prob = p1;
					}
					break;
				case IDLE: { 	// idle meaning going from state s to state s again
					switch (edge_position) // depending on which edge the state is (if any), the idle probability will be different
					{
						case RIGHT: prob = p1 + p2;	// if at the right edge, that means that agent is trying to go off the grid
								break;
						case DIAG_UR: prob = p1 + p2;
								break;
						case DIAG_DR: prob = p1 + p2;
								break;
						default: prob = p2;
								break;
					}
				}
				case UP: if (edge_position == RIGHT) prob = p3;
						break;
				case DOWN: if (edge_position == RIGHT) prob = p3;
						break;
				case DIAG_UR: prob = p3; 
						break;
				case DIAG_DR: prob = p3; 
						break;
			};
			break;
	};
	return prob;
}

void policy_evaluation() {
	double delta;
	int count = 0;
	do {
		vector<double> V_old (V);
		delta = 0.0;
		for (int s = 0; s < 16; s++) {  // for all s 
			double v = V_old[s];
			double sum1 = 0.0;
			for (int a = 0; a < 4; a++) { // for all actions
				double policy_a_given_s = policies[s][a];
				double sum2 = 0.0;
				for (int s_prime = 0; s_prime < 16; s_prime++) { // for all next states
					double prob = transition_probabilities[s][a][s_prime];
					double reward = rewards[a];
					sum2 = sum2 + (prob * (reward + gamma_const * V_old[s_prime]));
				}
				sum1 = sum1 + (policy_a_given_s * sum2);
			}
			
			V[s] = sum1;
			delta = max(delta, abs(v - V[s]));
		}
		count++;
	} while (delta > theta);
}

bool policy_improvement() {
	bool stable = true;
	
	for (int s = 0; s < 16; s++) {
		int max = 0;
		vector<double> action_val(4);
		for (int a = 0; a < 4; a++) {
			double sum_s_prime = 0.0;
			for (int s_prime = 0; s_prime < 16; s_prime++) {
				double prob = transition_probabilities[s][a][s_prime];
				double reward = rewards[a];
				sum_s_prime = sum_s_prime + (prob * (reward + gamma_const * V[s_prime]));
			}
			action_val[a] = sum_s_prime;
		}

		max = max_index(action_val); 

		if (!(abs(policies[s][max] - 1.0) <= 0.00001)) {
			stable = false;
			for (int i = 0; i < 4; i++) {
				if (i == max) {
					policies[s][i] = 1.0;
				} else {
					policies[s][i] = 0.0;
				}
			}
		}

	}
	return stable;
}

void policy_iteration() {
	bool policy_stable = false;
	int iterations = 0;
	initialize_probabilities(); // calculating transition probabilities
	do {
		clock_t start;
		start = clock();
		policy_evaluation();
		policy_stable = policy_improvement();
		iterations++;
		cout << "Iteration " << iterations << ", elapsed time: " << (clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << endl;
	} while (!policy_stable);
	cout << "Total iterations: " << iterations << endl;
}

int main(int argc, char *argv[])
{	
	if (argc > 2 && strcmp(argv[1], "-m") == 0) {
		p1 = atof(argv[2]);
		p2 = atof(argv[3]);

		for (int i = 0; i < 4; i++) {
			rewards[i] = atof(argv[4 + i]);
		}

	}

	policy_iteration();

	cout << "Optimal policy: " << endl;
	for (int s = 0; s < 16; s++) {
		char sign;
		switch (max_index(policies[s])) {
			case UP: sign = '^'; break;
			case DOWN: sign = 'v'; break;
			case LEFT: sign = '<'; break;
			case RIGHT: sign = '>'; break;
		}
		if ( s == 0 || s == 15) sign = ' ';
		printf("[%c]", sign);
		if (s == 3 || s == 7 || s == 11 || s == 15) { 
			cout << endl;
		}
	}

}

