CWD:=$(shell pwd)

SRC_DIR:=$(CWD)/src
TEST_SRC_DIR:=$(CWD)/src/test
BIN_DIR:=$(CWD)/bin
OBJ_DIR:=$(CWD)/obj
TEST_OBJ_DIR:=$(CWD)/obj/test

CXX:=g++
CXXFLAGS:=-std=c++11

INCLUDE_FLAGS:= -I$(SRC_DIR) -I$(TEST_SRC_DIR)

PROBABILITY_DEPS := $(SRC_DIR)/probability.hpp $(SRC_DIR)/reference.hpp $(SRC_DIR)/allele.hpp $(SRC_DIR)/input_haplotype.hpp $(SRC_DIR)/penalty_set.hpp $(SRC_DIR)/delay_multiplier.hpp $(SRC_DIR)/math.hpp $(SRC_DIR)/DP_map.hpp $(SRC_DIR)/row_set.hpp

CORE_OBJ := $(OBJ_DIR)/math.o $(OBJ_DIR)/reference.o $(OBJ_DIR)/probability.o $(OBJ_DIR)/input_haplotype.o $(OBJ_DIR)/delay_multiplier.o $(OBJ_DIR)/DP_map.o $(OBJ_DIR)/penalty_set.o $(OBJ_DIR)/allele.o $(OBJ_DIR)/row_set.o

TREE_OBJ := $(OBJ_DIR)/haplotype_state_node.o $(OBJ_DIR)/haplotype_state_tree.o $(OBJ_DIR)/haplotype_manager.o $(OBJ_DIR)/scored_node.o $(OBJ_DIR)/set_of_extensions.o $(OBJ_DIR)/reference_sequence.o

all : build_dirs speed speed_tree tests tree_tests interface

build_dirs:
	mkdir -p obj/test && mkdir -p bin

speed : $(TEST_OBJ_DIR)/speed.o $(CORE_OBJ)
	$(CXX) $(CXXFLAGS) $(TEST_OBJ_DIR)/speed.o $(CORE_OBJ) -o $(BIN_DIR)/speed

tests : $(TEST_OBJ_DIR)/test.o $(CORE_OBJ)
	$(CXX) $(CXXFLAGS) $(TEST_OBJ_DIR)/test.o $(CORE_OBJ) -o $(BIN_DIR)/tests

tree_tests : $(TEST_OBJ_DIR)/tree_tests.o $(CORE_OBJ) $(TREE_OBJ)
	$(CXX) $(CXXFLAGS) $(TEST_OBJ_DIR)/tree_tests.o $(CORE_OBJ) $(TREE_OBJ) -o $(BIN_DIR)/tree_tests

speed_tree : $(TEST_OBJ_DIR)/speed_tree.o $(CORE_OBJ) $(TREE_OBJ)
	$(CXX) $(CXXFLAGS) $(TEST_OBJ_DIR)/speed_tree.o $(CORE_OBJ) $(TREE_OBJ) -o $(BIN_DIR)/speed_tree

interface : $(OBJ_DIR)/linhapexample.o $(OBJ_DIR)/interface.o
	$(CXX) $(CXXFLAGS) $(OBJ_DIR)/linhapexample.o $(OBJ_DIR)/interface.o $(CORE_OBJ) $(TREE_OBJ) -o $(BIN_DIR)/linhapexample

clean:
	rm -f $(BIN_DIR)/* $(OBJ_DIR)/*.o $(TEST_OBJ_DIR)/*.o

$(OBJ_DIR)/allele.o : $(SRC_DIR)/allele.cpp $(SRC_DIR)/allele.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) -c $(SRC_DIR)/allele.cpp -o $(OBJ_DIR)/allele.o

$(OBJ_DIR)/linhapexample.o : $(SRC_DIR)/linhapexample.c $(SRC_DIR)/interface.h $(SRC_DIR)/haplotype_manager.hpp $(SRC_DIR)/reference_sequence.hpp $(SRC_DIR)/set_of_extensions.hpp $(SRC_DIR)/haplotype_state_tree.hpp $(SRC_DIR)/haplotype_state_node.hpp $(PROBABILITY_DEPS)
	gcc $(INCLUDE_FLAGS) -c $(SRC_DIR)/linhapexample.c -o $(OBJ_DIR)/linhapexample.o

$(OBJ_DIR)/interface.o : $(SRC_DIR)/interface.c $(SRC_DIR)/interface.h $(SRC_DIR)/haplotype_manager.hpp $(SRC_DIR)/reference_sequence.hpp $(SRC_DIR)/set_of_extensions.hpp $(SRC_DIR)/haplotype_state_tree.hpp $(SRC_DIR)/haplotype_state_node.hpp $(PROBABILITY_DEPS)
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) -c $(SRC_DIR)/interface.c -o $(OBJ_DIR)/interface.o

$(OBJ_DIR)/haplotype_state_node.o : $(SRC_DIR)/haplotype_state_node.cpp $(SRC_DIR)/haplotype_state_node.hpp $(PROBABILITY_DEPS)
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) -c $(SRC_DIR)/haplotype_state_node.cpp -o $(OBJ_DIR)/haplotype_state_node.o

$(OBJ_DIR)/haplotype_state_tree.o : $(SRC_DIR)/haplotype_state_tree.cpp $(SRC_DIR)/haplotype_state_tree.hpp $(SRC_DIR)/haplotype_state_node.hpp $(PROBABILITY_DEPS)
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) -c $(SRC_DIR)/haplotype_state_tree.cpp -o $(OBJ_DIR)/haplotype_state_tree.o

$(OBJ_DIR)/haplotype_manager.o : $(SRC_DIR)/haplotype_manager.cpp $(SRC_DIR)/haplotype_manager.hpp $(SRC_DIR)/reference_sequence.hpp $(SRC_DIR)/set_of_extensions.hpp $(SRC_DIR)/haplotype_state_tree.hpp $(SRC_DIR)/haplotype_state_node.hpp $(PROBABILITY_DEPS)
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) -c $(SRC_DIR)/haplotype_manager.cpp -o $(OBJ_DIR)/haplotype_manager.o

$(TEST_OBJ_DIR)/speed_tree.o : $(TEST_SRC_DIR)/speed_tree.cpp $(SRC_DIR)/haplotype_manager.hpp $(SRC_DIR)/reference_sequence.hpp $(SRC_DIR)/set_of_extensions.hpp $(SRC_DIR)/haplotype_state_tree.hpp $(SRC_DIR)/haplotype_state_node.hpp $(PROBABILITY_DEPS)
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) -c $(TEST_SRC_DIR)/speed_tree.cpp -o $(TEST_OBJ_DIR)/speed_tree.o

$(OBJ_DIR)/delay_multiplier.o : $(SRC_DIR)/delay_multiplier.cpp $(SRC_DIR)/delay_multiplier.hpp $(SRC_DIR)/math.hpp $(SRC_DIR)/DP_map.hpp $(SRC_DIR)/row_set.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) -c $(SRC_DIR)/delay_multiplier.cpp -o $(OBJ_DIR)/delay_multiplier.o

$(OBJ_DIR)/DP_map.o : $(SRC_DIR)/DP_map.cpp $(SRC_DIR)/math.hpp $(SRC_DIR)/DP_map.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) -c $(SRC_DIR)/DP_map.cpp -o $(OBJ_DIR)/DP_map.o

$(OBJ_DIR)/input_haplotype.o : $(SRC_DIR)/input_haplotype.cpp $(SRC_DIR)/input_haplotype.hpp $(SRC_DIR)/reference.hpp $(SRC_DIR)/allele.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) -c $(SRC_DIR)/input_haplotype.cpp -o $(OBJ_DIR)/input_haplotype.o

$(OBJ_DIR)/math.o : $(SRC_DIR)/math.cpp $(SRC_DIR)/math.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) -c $(SRC_DIR)/math.cpp -o $(OBJ_DIR)/math.o

$(OBJ_DIR)/probability.o : $(SRC_DIR)/probability.cpp $(PROBABILITY_DEPS)
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) -c $(SRC_DIR)/probability.cpp -o $(OBJ_DIR)/probability.o

$(OBJ_DIR)/penalty_set.o : $(SRC_DIR)/penalty_set.cpp $(SRC_DIR)/penalty_set.hpp $(SRC_DIR)/math.hpp $(SRC_DIR)/DP_map.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) -c $(SRC_DIR)/penalty_set.cpp -o $(OBJ_DIR)/penalty_set.o

$(OBJ_DIR)/reference.o : $(SRC_DIR)/reference.cpp $(SRC_DIR)/reference.hpp $(SRC_DIR)/allele.hpp $(SRC_DIR)/row_set.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) -c $(SRC_DIR)/reference.cpp -o $(OBJ_DIR)/reference.o

$(OBJ_DIR)/reference_sequence.o : $(SRC_DIR)/reference_sequence.cpp $(SRC_DIR)/reference_sequence.hpp $(SRC_DIR)/reference.hpp $(SRC_DIR)/allele.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) -c $(SRC_DIR)/reference_sequence.cpp -o $(OBJ_DIR)/reference_sequence.o

$(OBJ_DIR)/row_set.o : $(SRC_DIR)/row_set.cpp $(SRC_DIR)/row_set.hpp $(SRC_DIR)/allele.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) -c $(SRC_DIR)/row_set.cpp -o $(OBJ_DIR)/row_set.o

$(OBJ_DIR)/scored_node.o : $(SRC_DIR)/scored_node.cpp $(SRC_DIR)/scored_node.hpp $(SRC_DIR)/haplotype_state_node.hpp $(PROBABILITY_DEPS)
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) -c $(SRC_DIR)/scored_node.cpp -o $(OBJ_DIR)/scored_node.o

$(OBJ_DIR)/set_of_extensions.o : $(SRC_DIR)/set_of_extensions.cpp $(SRC_DIR)/set_of_extensions.hpp  $(PROBABILITY_DEPS)
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) -c $(SRC_DIR)/set_of_extensions.cpp -o $(OBJ_DIR)/set_of_extensions.o

$(TEST_OBJ_DIR)/test.o : $(TEST_SRC_DIR)/test.cpp $(PROBABILITY_DEPS)
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) -c $(TEST_SRC_DIR)/test.cpp -o $(TEST_OBJ_DIR)/test.o

$(TEST_OBJ_DIR)/tree_tests.o : $(TEST_SRC_DIR)/tree_tests.cpp $(SRC_DIR)/haplotype_manager.hpp $(SRC_DIR)/reference_sequence.hpp $(SRC_DIR)/set_of_extensions.hpp $(SRC_DIR)/haplotype_state_tree.hpp $(SRC_DIR)/haplotype_state_node.hpp $(PROBABILITY_DEPS)
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) -c $(TEST_SRC_DIR)/tree_tests.cpp -o $(TEST_OBJ_DIR)/tree_tests.o

$(TEST_OBJ_DIR)/speed.o : $(TEST_SRC_DIR)/speed.cpp $(PROBABILITY_DEPS)
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) -c $(TEST_SRC_DIR)/speed.cpp -o $(TEST_OBJ_DIR)/speed.o
