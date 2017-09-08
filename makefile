SRC_DIR:=src
TEST_SRC_DIR:=src/test
BIN_DIR:=bin
OBJ_DIR:=obj
TEST_OBJ_DIR:=obj/test

CXX:=g++
CXXFLAGS:=-std=c++11

CWD:=$(shell pwd)

PROBABILITY_DEPS := $(SRC_DIR)/probability.hpp $(SRC_DIR)/reference.hpp $(SRC_DIR)/allele.hpp $(SRC_DIR)/input_haplotype.hpp $(SRC_DIR)/penalty_set.hpp $(SRC_DIR)/delay_multiplier.hpp $(SRC_DIR)/math.hpp $(SRC_DIR)/DP_map.hpp

VCF_DEPS := -I deps/vcflib/include -I deps/vcflib/tabixpp/htslib -Ldeps/vcflib/lib -Ldeps/vcflib/tabixpp/htslib -lvcflib -lhts -lz -lpthread -lm

CORE_OBJ := math.o reference.o probability.o input_haplotype.o delay_multiplier.o DP_map.o penalty_set.o allele.o
TREE_OBJ := haplotype_state_node.o haplotype_state_tree.o haplotype_manager.o vcf_manager.o scored_node.o set_of_extensions.o reference_sequence.o

speed : speed.o $(CORE_OBJ)
	$(CXX) $(CXXFLAGS) speed.o $(CORE_OBJ) -o speed

tests : test.o $(CORE_OBJ)
	$(CXX) $(CXXFLAGS) test.o $(CORE_OBJ) -o tests

tree_tests: tree_tests.o $(CORE_OBJ) $(TREE_OBJ)
	$(CXX) $(CXXFLAGS) tree_tests.o $(CORE_OBJ) $(TREE_OBJ) -o tree_tests

$(OBJ_DIR)/allele.o : $(SRC_DIR)/allele.cpp $(SRC_DIR)/allele.hpp
	$(CXX) $(CXXFLAGS) -c allele.cpp

$(OBJ_DIR)/haplotype_manager.o : $(SRC_DIR)/haplotype_manager.cpp $(SRC_DIR)/haplotype_manager.hpp
	$(CXX) $(CXXFLAGS) -c haplotype_manager.cpp

$(OBJ_DIR)/haplotype_state_node.o : $(SRC_DIR)/haplotype_state_node.cpp $(SRC_DIR)/haplotype_state_node.hpp $(PROBABILITY_DEPS)
	$(CXX) $(CXXFLAGS) -c haplotype_state_node.cpp

$(OBJ_DIR)/haplotype_state_tree.o : $(SRC_DIR)/haplotype_state_tree.cpp $(SRC_DIR)/haplotype_state_tree.hpp $(SRC_DIR)/haplotype_state_node.hpp $(PROBABILITY_DEPS)
	$(CXX) $(CXXFLAGS) -c haplotype_state_tree.cpp

$(OBJ_DIR)/haplotype_manager.o : $(SRC_DIR)/haplotype_manager.cpp $(SRC_DIR)/haplotype_manager.hpp $(SRC_DIR)/reference_sequence.hpp $(SRC_DIR)/set_of_extensions.hpp $(SRC_DIR)/haplotype_state_tree.hpp $(SRC_DIR)/haplotype_state_node.hpp $(PROBABILITY_DEPS)
	$(CXX) $(CXXFLAGS) -c haplotype_manager.cpp

$(OBJ_DIR)/delay_multiplier.o : $(SRC_DIR)/delay_multiplier.cpp $(SRC_DIR)/delay_multiplier.hpp $(SRC_DIR)/math.hpp $(SRC_DIR)/DP_map.hpp
	$(CXX) $(CXXFLAGS) -c delay_multiplier.cpp

$(OBJ_DIR)/DP_map.o : $(SRC_DIR)/DP_map.cpp $(SRC_DIR)/math.hpp $(SRC_DIR)/DP_map.hpp
	$(CXX) $(CXXFLAGS) -c DP_map.cpp

$(OBJ_DIR)/input_haplotype.o : $(SRC_DIR)/input_haplotype.cpp $(SRC_DIR)/input_haplotype.hpp $(SRC_DIR)/reference.hpp $(SRC_DIR)/allele.hpp
	$(CXX) $(CXXFLAGS) -c input_haplotype.cpp

$(OBJ_DIR)/math.o : $(SRC_DIR)/math.cpp $(SRC_DIR)/math.hpp
	$(CXX) $(CXXFLAGS) -c math.cpp

$(OBJ_DIR)/probability.o : $(SRC_DIR)/probability.cpp $(PROBABILITY_DEPS)
	$(CXX) $(CXXFLAGS) -c probability.cpp

$(OBJ_DIR)/penalty_set.o : $(SRC_DIR)/penalty_set.cpp $(SRC_DIR)/penalty_set.hpp $(SRC_DIR)/math.hpp $(SRC_DIR)/DP_map.hpp
	$(CXX) $(CXXFLAGS) -c penalty_set.cpp

$(OBJ_DIR)/reference.o : $(SRC_DIR)/reference.cpp $(SRC_DIR)/reference.hpp $(SRC_DIR)/allele.hpp
	$(CXX) $(CXXFLAGS) -c reference.cpp

$(OBJ_DIR)/reference_sequence.o : $(SRC_DIR)/reference_sequence.cpp $(SRC_DIR)/reference_sequence.hpp $(SRC_DIR)/reference.hpp $(SRC_DIR)/allele.hpp
	$(CXX) $(CXXFLAGS) -c reference_sequence.cpp

$(OBJ_DIR)/scored_node.o : $(SRC_DIR)/scored_node.cpp $(SRC_DIR)/scored_node.hpp $(SRC_DIR)/haplotype_state_node.hpp $(PROBABILITY_DEPS)
	$(CXX) $(CXXFLAGS) -c scored_node.cpp

$(OBJ_DIR)/set_of_extensions.o : $(SRC_DIR)/set_of_extensions.cpp $(SRC_DIR)/set_of_extensions.hpp  $(PROBABILITY_DEPS)
	$(CXX) $(CXXFLAGS) -c set_of_extensions.cpp

$(OBJ_DIR)/vcf_manager.o : $(SRC_DIR)/vcf_manager.cpp $(SRC_DIR)/vcf_manager.hpp $(SRC_DIR)/reference.hpp $(SRC_DIR)/allele.hpp
	$(CXX) $(CXXFLAGS) $(VCF_DEPS) -c vcf_manager.cpp

$(TEST_OBJ_DIR)/test.o : $(TEST_SRC_DIR)/test.cpp

$(TEST_OBJ_DIR)/tree_tests.o : $(TEST_SRC_DIR)/tree_tests.cpp
