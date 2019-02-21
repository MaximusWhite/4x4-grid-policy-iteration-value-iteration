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
vector<double> rewards = { -1.0, -1.0, -1.0, -1.0};

vector<double> V(16, 0.0); // initial V(s) = 0, for all s
vector<vector<vector<double>>> transition_probabilities(16, vector<vector<double>>(4, vector<double>(16, 0.0))); // s, a, s'
vector<vector<double>> policies(16, vector<double>(4, 0.25)); // policy is equiprobable for each action in any given state 

int get_direction(int state, int state_prime);
void initialize_probabilities();
int is_on_edge(int state);
double action_probability(int state, int state_prime, int action);
int max_index(vector<double> arr);
void value_iteration();


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

void initialize_V() {
    for (int i = 1; i < V.size() - 1; i++) {
        V[i] = rand();
    }
}

void value_iteration() {
    double delta;
	int count = 0;
    vector<int> best_choises(16, 0);
	do {
		delta = 0.0;
		for (int s = 0; s < 16; s++) {  // for all s 
			double v = V[s];
			double sum1 = 0.0;
            vector<double> sums(4, -5);
			for (int a = 0; a < 4; a++) { // for all actions
				double policy_a_given_s = policies[s][a];
				double sum2 = 0.0;
                double reward = rewards[a];
				for (int s_prime = 0; s_prime < 16; s_prime++) { // for all next states
					double prob = transition_probabilities[s][a][s_prime];
					sum2 = sum2 + (prob * (reward + gamma_const * V[s_prime]));
                    printf("S = %d, a = %d, s' = %d, sum2 = %f, prob = %f\n", s, a, s_prime, sum2, prob);

				}
				sum1 = sum1 + (policy_a_given_s * sum2);

                // printf("S = %d, a = %d, sum2 = %f, sum1 = %f\n", s, a, sum2, sum1);
                sums[a] = sum1;
			}
    
            int max_a = max_index(sums);
            best_choises[s] = max_a;
			V[s] = sums[max_a];
			delta = max(delta, abs(v - V[s]));
		}
	} while (delta > theta);

    for (int s = 0; s < 16; s++) {      // argmax for policy per state 
        for (int a = 0; a < 4; a++) {
            if (best_choises[s] == a) {
                policies[s][a] = 1.0;
            } else {
                policies[s][a] = 0.0;
            }
        }
    }

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
    initialize_probabilities();
    value_iteration();

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

