#define CATCH_CONFIG_MAIN

#include <cmath> 
#include "haplotype_state_tree.hpp"
#include "haplotype_state_node.hpp"
#include "probability.hpp"
#include "reference_sequence.hpp"
#include "haplotype_manager.hpp"
#include "catch.hpp"

using namespace std;

linearReferenceStructure build_ref(string ref_seq, vector<size_t> positions) {
  vector<alleleValue> ref_values;
  for(size_t i = 0; i < positions.size(); i++) {
    ref_values.push_back(char_to_allele(ref_seq[positions[i]]));
  }
  return linearReferenceStructure(positions, ref_seq.length(), ref_values);
}

TEST_CASE( "Node construction and destruction works as intended", "[node]") {
  haplotypeStateNode* n1 = new haplotypeStateNode();
  REQUIRE(n1->is_leaf());
  REQUIRE(n1->is_abandoned_stem());
  haplotypeStateNode* n2 = n1->add_child(A);
  haplotypeStateNode* n3 = n1->add_child(C);
  REQUIRE(n1->number_of_children() == 2);
  REQUIRE(n2->get_allele() == A);
  REQUIRE(n3->get_allele() == C);
  REQUIRE(n2 == n1->get_child(A));
  REQUIRE(n2->get_parent() == n1);
  n1->remove_child(A);
  REQUIRE(n1->number_of_children() == 1);
  n1->remove_child(C);
  REQUIRE(n1->number_of_children() == 0);
}

/*
TEST_CASE( "Tree construction works as intended", "[tree]" ) {
  
}
*/

TEST_CASE( "Pruning functions work", "[tree]" ) {
  linearReferenceStructure reference = build_ref("AAAA", {0,1,2,3});
  vector<vector<alleleValue> > haplotypes = {
    {A, A, A, A},
    {A, T, A, A},
    {A, A, T, A},
    {A, A, A, T},
    {A, A, A, A}
  };
  haplotypeCohort cohort = haplotypeCohort(haplotypes, &reference);
  penaltySet penalties = penaltySet(-6, -9, 3);
  
  string ref_string = "AAAA";
  string read_string = "AAAA";
  vector<size_t> read_sites = {0,1,2,3};
  
  haplotypeStateTree hsTree = haplotypeStateTree(&reference, &penalties, 
              &cohort);
  haplotypeStateNode* n;
  n = hsTree.root->add_child(A);
  n = n->add_child(C);
  n = n->add_child(T);
  n = n->add_child(G);
  REQUIRE(hsTree.root->number_of_children() == 1);
  hsTree.remove_node_and_unshared_ancestors(n);
  REQUIRE(hsTree.root->number_of_children() == 0);
  haplotypeManager hap_manager = haplotypeManager(
            &reference,
            &cohort,
            &penalties,
            ref_string.c_str(),
            read_sites,
            read_string.c_str(),
            0);
}


TEST_CASE( "Tree navigation works as intended", "[tree]" ) {
  linearReferenceStructure reference = build_ref("AAA", {0,1,2});
  vector<vector<alleleValue> > haplotypes = {
    {A, A, A},
    {C, C, C},
    {T, T, T},
    {G, G, G},
    {gap, gap, gap}
  };
  haplotypeCohort cohort = haplotypeCohort(haplotypes, &reference);
  penaltySet penalties = penaltySet(-6, -9, 3);
  
  haplotypeStateTree hsTree = haplotypeStateTree(&reference, &penalties, 
              &cohort);
  haplotypeStateNode* n;
  n = hsTree.root->add_child(A);
  n = n->add_child(C);
  n = n->add_child(T);
  
  // Can retrieve path from node
  vector<alleleValue> test = hsTree.state_to_alleles(n);
  REQUIRE(test.size() == 3);
  REQUIRE(test[0] == A);
  REQUIRE(test[1] == C);
  REQUIRE(test[2] == T);
  
  n = hsTree.root->add_child(T);
  n = n->add_child(G);
  n = n->add_child(A);
    
  // Can retrieve node from path
  vector<alleleValue> query = {T, G, A};
  haplotypeStateNode* test2 = hsTree.alleles_to_state(query);
  REQUIRE(test2 == n);
}


TEST_CASE( "Tree interfaces with probability DP state matrices", "[tree][DP]") {
  linearReferenceStructure reference = build_ref("AAAA", {0,1,2,3});
  vector<vector<alleleValue> > haplotypes = {
    {A, A, A, A},
    {A, T, A, A},
    {A, A, T, A},
    {A, A, A, T},
    {A, A, A, A}
  };
  haplotypeCohort cohort = haplotypeCohort(haplotypes, &reference);
  penaltySet penalties = penaltySet(-6, -9, 3);
  string ref_string = "AAAA";
  string read_string = "AAAA";
  vector<size_t> read_sites = {0,1,2,3};
    
  SECTION( "Copying preserves correctly-extending states" ) {
    haplotypeStateTree hsTree = haplotypeStateTree(&reference, &penalties, 
                &cohort);
    haplotypeStateNode* n;
    haplotypeMatrix correct_state = haplotypeMatrix(
                &reference,
                &penalties,
                &cohort);
    correct_state.initialize_probability_at_site(0, A);
    n = hsTree.root->add_child(A);
    n->state->initialize_probability_at_site(0, A);
    REQUIRE(correct_state.S == n->state->S);
    n = n->add_child_copying_state(C);
    correct_state.extend_probability_at_site(1, C);
    n->state->extend_probability_at_site(1, C);
    REQUIRE(correct_state.S == n->state->S);
  }
  
  SECTION( "HaplotypeManager maintains correct states for simple sets of sites" ) {
    haplotypeMatrix correct_state = haplotypeMatrix(
                &reference,
                &penalties,
                &cohort);
    
    haplotypeManager hap_manager = haplotypeManager(
              &reference,
              &cohort,
              &penalties,
              ref_string.c_str(),
              read_sites,
              read_string.c_str(),
              0);

    hap_manager.build_entire_tree(0);

    haplotypeStateNode* n = hap_manager.get_tree()->root;
    correct_state.initialize_probability_at_site(0, A);
    n = n->get_child(A);
    REQUIRE(correct_state.S == n->state->S);
    correct_state.extend_probability_at_site(1, C);
    n = n->get_child(C);
    REQUIRE(correct_state.S == n->state->S);
    correct_state.extend_probability_at_site(2, T);
    n = n->get_child(T);
    REQUIRE(correct_state.S == n->state->S);
    correct_state.extend_probability_at_site(3, G);
    n = n->get_child(G);
    REQUIRE(correct_state.S == n->state->S);
  }
}

TEST_CASE( "Haplotype manager performs correct calculations in the presence of unshared sites", "[manager]" ) {
  string ref_refstring = "AAAAAAAAAAAAAAAAAAAA";
  string read_refstring = "AAAAAAAAAAAAAAAAAA";
  vector<size_t> ref_sites = {1,3,4,8,13,19};
  vector<size_t> read_sites = {3,5,6,13,17};
  linearReferenceStructure reference = build_ref(ref_refstring, ref_sites);
  vector<vector<alleleValue> > haplotypes = {
   {A, A, A, A, A, A},
   {A, T, A, A, A, A},
   {A, A, T, A, A, A},
   {A, A, A, T, A, A},
   {A, A, A, A, T, A},
   {A, A, A, A, A, T}
  };
  haplotypeCohort cohort = haplotypeCohort(haplotypes, &reference);
  penaltySet penalties = penaltySet(-6, -9, 3);
  
  SECTION( "Haplotype manager performs correct calculations without pruning") {
    haplotypeManager hap_manager = haplotypeManager(
                &reference, 
                &cohort, 
                &penalties, 
                ref_refstring.c_str(),
                read_sites, 
                read_refstring.c_str(),
                2);
    
  }
  
  SECTION( "Haplotype manager performs correct calculations with pruning") {
    haplotypeManager hap_manager = haplotypeManager(
                &reference, 
                &cohort, 
                &penalties, 
                ref_refstring.c_str(),
                read_sites, 
                read_refstring.c_str(),
                2);
    
  }
}

TEST_CASE( "Haplotype Manager can translate positions and identify and translate sites" , "[manager]") {
  string ref_refstring = "AAAAAAAAAAAAAAAAAAAA";
  string read_refstring = "AAAAAAAAAAAAAAAAAA";
  vector<size_t> ref_sites = {1,3,4,8,13,19};
  vector<size_t> read_sites = {3,5,6,13,17};
  linearReferenceStructure reference = build_ref(ref_refstring, ref_sites);
  vector<vector<alleleValue> > haplotypes = {
   {A, A, A, A, A, A}
  };
  haplotypeCohort cohort = haplotypeCohort(haplotypes, &reference);
  penaltySet penalties = penaltySet(-6, -9, 3);
  referenceSequence ref_refseq = referenceSequence(ref_refstring);
  
  haplotypeManager hap_manager = haplotypeManager(
              &reference, 
              &cohort, 
              &penalties, 
              ref_refstring.c_str(),
              read_sites, 
              read_refstring.c_str(),
              2);
  
  REQUIRE(hap_manager.length() == 18);
  REQUIRE(hap_manager.read_sites() == 5);
  REQUIRE(hap_manager.shared_sites() == 2);
  REQUIRE(hap_manager.ref_position(4) == 6);
  REQUIRE(hap_manager.read_position(3) == 1);
  REQUIRE(hap_manager.get_read_site_read_position(1) == 5);
  REQUIRE(hap_manager.get_read_site_ref_position(1) == 7);
  REQUIRE(hap_manager.index_among_shared_sites(4) == 1);
  REQUIRE(hap_manager.index_among_read_only_sites(3) == 2);
  REQUIRE(hap_manager.get_shared_site_read_index(1) == 4);
  REQUIRE(hap_manager.get_shared_site_ref_index(1) == 5);
  REQUIRE(hap_manager.get_ref_site_below_read_site(2) == 3);
  REQUIRE(hap_manager.get_ref_site_below_read_site(3) == 4);
  
  //TODO test penalties
}
