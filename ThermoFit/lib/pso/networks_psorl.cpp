#pragma once
#include "networks_psorl.hpp"

Network::Network() {
	resetNetwork();
}

void Network::resetNetwork() {
	A1.resize(64, 0.0);
	A2.resize(64, 0.0);
	A3.resize(64, 0.0);
	A4.resize(25, 0.0);
	A5.resize(15, 0.0);
}

void Network::leakyRELU(double* x) {
	if (*x < 0.0) {
		*x = 0.01 * *x;
	}
}

void Network::tanh(double* x) {
	*x = std::tanh(*x);
}

std::vector<double> Network::forward(std::vector<double> input) {
	
	for (int i = 0; i < A1.size(); i++) {
		A1[i] = 0.0;
		for (int j = 0; j < input.size(); j++) {
			A1[i] += net_W1[i][j] * input[j];
		}
		A1[i] += net_B1[i][0];
		leakyRELU(&A1[i]);
	}

	for (int i = 0; i < A2.size(); i++) {
		A2[i] = 0.0;
		for (int j = 0; j < A1.size(); j++) {
			A2[i] += net_W2[i][j] * A1[j];
		}
		A2[i] += net_B2[i][0];
		leakyRELU(&A2[i]);
	}

	for (int i = 0; i < A3.size(); i++) {
		A3[i] = 0.0;
		for (int j = 0; j < A2.size(); j++) {
			A3[i] += net_W3[i][j] * A2[j];
		}
		A3[i] += net_B3[i][0];
		leakyRELU(&A3[i]);
	}

	for (int i = 0; i < A4.size(); i++) {
		A4[i] = 0.0;
		for (int j = 0; j < A3.size(); j++) {
			A4[i] += net_W4[i][j] * A3[j];
		}
		A4[i] += net_B4[i][0];
		leakyRELU(&A4[i]);
	}

	for (int i = 0; i < A5.size(); i++) {
		A5[i] = 0.0;
		for (int j = 0; j < A4.size(); j++) {
			A5[i] += net_W5[i][j] * A4[j];
		}
		A5[i] += net_B5[i][0];
		tanh(&A5[i]);
	}
	
	return A5;
}