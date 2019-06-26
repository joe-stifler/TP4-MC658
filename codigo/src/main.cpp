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

#include <DCMSTP.h>
#include <DCMSTPLagrangean.h>
#include <DCMSTPMetaheuristic.h>

int main(int argc, const char **argv) {
	if (argc == 4) {
		int d;
		int n, m;
		int u, v, w;
		std::string line;
		std::stringstream input;
		std::ifstream instanceFile;
		clock_t initialTime = clock();
		std::unique_ptr<DCMSTP> solver;
		int maximumTime = atoi(argv[2]);

		instanceFile.open(argv[1], std::ios::in | std::ios::binary);

		if (!instanceFile.is_open()) {
    		printf("ERROR: Not possible to open instance file.\n");
			return 1;
		}

		/* Reads graph */
		std::getline(instanceFile, line);

		input = std::stringstream(line);
		input >> n >> m;

		if (argv[3][0] == 'l') {
			/* Lagrangian heuristic */
			solver.reset(new DCMSTPLagrangean(n, maximumTime, initialTime));
		} else if (argv[3][0] == 'm') {
			/* Metaheuristic */
			solver.reset(new DCMSTPMetaheuristic(n, maximumTime, initialTime));
		} else {
			printf("ERROR: You should correctly especify the fourth param, where it "
		           "must be either 'l' for lagrangean heuristic or 'm' for metaheuristic.\n");

			return 1;
		}

		if (solver.get()) {
			/* Reads all edges */
			while (m--) {
				std::getline(instanceFile, line);

                input = std::stringstream(line);
				input >> u >> v >> w;

				solver->addEdge(u - 1, v - 1, w);
			}

			/* Reads vertex degree constraints */
			while (n--) {
				std::getline(instanceFile, line);

                input = std::stringstream(line);
				input >> u >> d;

				solver->setVertexMaxDegree(u - 1, d);
			}

			solver->solve();

			solver->saveBestEdges(std::string(argv[1]));

            solver->printSolution(std::string(argv[1]));
		} else {
			printf("Some error occurred during DCMSTP memory class allocation.\n");

			return 1;
		}

		return 0;
	}

	printf("ERROR: you should pass 4 parameters: executable, instance file, "
		   "time param and method param. Instead, %d was/were passed!\n", argc);

	return 1;
}
