######################################################################################################################################
# HEADER - PRACTICAL WORK 4
######################################################################################################################################
#
# Discipline: MC658
# PED: Natanael Ramos
# Professor: Cid C. de Souza
# Data of creation: May 31, 2019
# Author (RA 176665): Jose Ribeiro Neto <j176665@dac.unicamp.br>
# Author (RA 171119): Felipe Lopes De Mello <f171119@dac.unicamp.br>
#
######################################################################################################################################
# BASIC SETTINGS (objects, flags, etc.)
######################################################################################################################################

OBJ_DIR = obj
SRC_DIR = src

OBJS = $(OBJ_DIR)/main.o $(OBJ_DIR)/DCMSTP.o $(OBJ_DIR)/DCMSTPLagrangean.o $(OBJ_DIR)/DisjointSets.o $(OBJ_DIR)/DCMSTPMetaheuristic.o

BINARY_NAME = dcmstp-solver
LINKING_FLAGS = -O3 -std=c++11
COMPILATION_FLAGS = -c -O3 -std=c++11 -Iinclude

######################################################################################################################################
# COMPILATION
######################################################################################################################################
all: create_obj_dir $(OBJS)
		g++ $(OBJS) -o $(BINARY_NAME) $(LINKING_FLAGS)

create_obj_dir:
		mkdir -p $(OBJ_DIR)

$(OBJS): $(OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp
		g++ -c $< -o $@ $(COMPILATION_FLAGS)

######################################################################################################################################
clean:
		rm -rf $(OBJ_DIR) $(BINARY_NAME)
