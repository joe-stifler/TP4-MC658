/**
 * File: main.cpp
 *
 * Discipline: MC658
 * PED: Natanael Ramos
 * Professor: Cid C. de Souza
 * Data of creation: May 31, 2019
 * Author (RA 176665): Jose Ribeiro Neto <j176665@dac.unicamp.br>
 * Author (RA 171119): Felipe Lopes De Mello <f171119@dac.unicamp.br>
 *
 **/

#include <cstdio>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iostream>

int main(int argc, const char **argv) {
	if (argc == 4) {
		clock_t initialTime = clock();

		long int maximumTime = atol(argv[2]);

		if (argv[3] == "l") {
			/* lagrangian heuristic */
		} else if (argv[3] == "m") {
			/* metaheuristic */
		} else {
			printf("ERROR: You should correctly especify the fourth param, where it must be either 'l' for lagrangan heuristic or 'm' for metaheuristic.\n");

			return 1;
		}

		std::ifstream instanceFile;
		instanceFile.open(argv[1], std::ios::in | std::ios::binary);

		if (!instanceFile.is_open()) {
    		printf("ERROR: Not possible to open instance file.\n");
			return 1;
		}

		std::string line;
		std::getline(instanceFile, line); /* Reads the params data (maximum number of nodes to explore). */

		int variable;
		std::stringstream test(line);
		test >> variable;

		return 0;
	}

	printf("ERROR: you should pass 4 parameters: executable, instance file, time param and method param. Instead, %d was/were passed!\n", argc);

	return 1;
}
