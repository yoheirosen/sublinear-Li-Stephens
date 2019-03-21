CWD:=$(shell pwd)

SRC_DIR:=$(CWD)/src
TEST_SRC_DIR:=$(CWD)/src/test
BIN_DIR:=$(CWD)/bin
OBJ_DIR:=$(CWD)/obj
TEST_OBJ_DIR:=$(CWD)/obj/test
LIB_DIR:= $(CWD)/lib
DEP_DIR:= $(CWD)/deps

CXX:=g++
CXXFLAGS:=-std=c++11

INCLUDE_FLAGS:= -I$(SRC_DIR) -I$(TEST_SRC_DIR) -I$(DEP_DIR)/htslib
LIBS := -L. -L$(DEP_DIR)/htslib/ -lhts -llzma -lbz2 -lz -lm -lpthread

LIBHTS := $(DEP_DIR)/htslib/libhts.a

PROBABILITY_DEPS := $(SRC_DIR)/probability.hpp $(SRC_DIR)/reference.hpp $(SRC_DIR)/allele.hpp $(SRC_DIR)/input_haplotype.hpp $(SRC_DIR)/penalty_set.hpp $(SRC_DIR)/delay_multiplier.hpp $(SRC_DIR)/math.hpp $(SRC_DIR)/DP_map.hpp $(SRC_DIR)/row_set.hpp

CORE_OBJ := $(OBJ_DIR)/reference.o $(OBJ_DIR)/probability.o $(OBJ_DIR)/input_haplotype.o $(OBJ_DIR)/delay_multiplier.o $(OBJ_DIR)/DP_map.o $(OBJ_DIR)/penalty_set.o $(OBJ_DIR)/allele.o $(OBJ_DIR)/row_set.o $(LIBHTS)

all : build_dirs tests libs serializer

build_dirs:
	if [ ! -d $(OBJ_DIR) ]; then mkdir -p $(OBJ_DIR); fi
	if [ ! -d $(TEST_OBJ_DIR) ]; then mkdir -p $(TEST_OBJ_DIR); fi
	if [ ! -d $(BIN_DIR) ]; then mkdir -p $(BIN_DIR); fi
	if [ ! -d $(LIB_DIR) ]; then mkdir -p $(LIB_DIR); fi
	if [ ! -d $(DEP_DIR) ]; then mkdir -p $(DEP_DIR); fi

tests : $(TEST_OBJ_DIR)/test.o $(CORE_OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $(BIN_DIR)/tests $(LIBS)

serializer : $(OBJ_DIR)/serialize_index.o $(CORE_OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $(BIN_DIR)/serializer $(LIBS)

profiler : $(OBJ_DIR)/profiler.o $(CORE_OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $(BIN_DIR)/profiler $(LIBS)

libs : $(LIB_DIR)/libsublinearLS.a $(CORE_OBJ)

clean:
	rm -f $(BIN_DIR)/* $(OBJ_DIR)/*.o $(TEST_OBJ_DIR)/*.o $(LIB_DIR)/*

$(LIB_DIR)/libsublinearLS.a : $(CORE_OBJ)
	ar rc $@ $^
	ranlib $@

$(OBJ_DIR)/allele.o : $(SRC_DIR)/allele.cpp $(SRC_DIR)/allele.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) $(LIBS) -c $< -o $@

$(OBJ_DIR)/delay_multiplier.o : $(SRC_DIR)/delay_multiplier.cpp $(SRC_DIR)/delay_multiplier.hpp $(SRC_DIR)/DP_map.hpp $(SRC_DIR)/row_set.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) $(LIBS) -c $< -o $@

$(OBJ_DIR)/DP_map.o : $(SRC_DIR)/DP_map.cpp $(SRC_DIR)/math.hpp $(SRC_DIR)/DP_map.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) $(LIBS) -c $< -o $@

$(OBJ_DIR)/input_haplotype.o : $(SRC_DIR)/input_haplotype.cpp $(SRC_DIR)/input_haplotype.hpp $(SRC_DIR)/reference.hpp $(SRC_DIR)/allele.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) $(LIBS) -c $< -o $@

$(OBJ_DIR)/probability.o : $(SRC_DIR)/probability.cpp $(PROBABILITY_DEPS)
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) $(LIBS) -c $< -o $@

$(OBJ_DIR)/penalty_set.o : $(SRC_DIR)/penalty_set.cpp $(SRC_DIR)/penalty_set.hpp $(SRC_DIR)/math.hpp $(SRC_DIR)/DP_map.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) $(LIBS) -c $< -o $@

$(OBJ_DIR)/reference.o : $(SRC_DIR)/reference.cpp $(SRC_DIR)/reference.hpp $(SRC_DIR)/allele.hpp $(SRC_DIR)/row_set.hpp $(LIBHTS)
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) $(LIBS) -c $< -o $@

$(OBJ_DIR)/reference_sequence.o : $(SRC_DIR)/reference_sequence.cpp $(SRC_DIR)/reference_sequence.hpp $(SRC_DIR)/reference.hpp $(SRC_DIR)/allele.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) $(LIBS) -c $< -o $@

$(OBJ_DIR)/row_set.o : $(SRC_DIR)/row_set.cpp $(SRC_DIR)/row_set.hpp $(SRC_DIR)/allele.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) $(LIBS) -c $< -o $@
	
$(TEST_OBJ_DIR)/test.o : $(TEST_SRC_DIR)/test.cpp $(PROBABILITY_DEPS)
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) $(LIBS) -c $< -o $@

$(OBJ_DIR)/serialize_index.o : $(SRC_DIR)/serialize_index.cpp $(SRC_DIR)/reference.hpp $(SRC_DIR)/allele.hpp $(SRC_DIR)/row_set.hpp $(LIBHTS)
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) $(LIBS) -c $< -o $@

$(OBJ_DIR)/profiler.o : $(SRC_DIR)/test/profiler.cpp $(PROBABILITY_DEPS)
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) $(LIBS) -c $< -o $@

$(LIBHTS) :
	cd deps/htslib && make
