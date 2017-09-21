#define CATCH_CONFIG_MAIN

#include <cmath> 
#include "haplotype_state_tree.hpp"
#include "haplotype_state_node.hpp"
#include "probability.hpp"
#include "reference_sequence.hpp"
#include "haplotype_manager.hpp"
#include "catch.hpp"
#include "vcf_manager.hpp"
#include <iostream>

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

/*
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
*/

TEST_CASE( "Tree navigation works as intended", "[tree][navigation]" ) {
  linearReferenceStructure reference = build_ref("AAA", {0,1,2});
  vector<vector<alleleValue> > haplotypes = {
    {A, A, A},
    {C, C, C},
    {T, T, T},
    {G, G, G},
    {gap, gap, gap}
  };
  haplotypeCohort cohort = haplotypeCohort(haplotypes, &reference);
  penaltySet penalties = penaltySet(-6, -9, haplotypes.size());
  
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
  penaltySet penalties = penaltySet(-6, -9, haplotypes.size());
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
    n->state = new haplotypeMatrix(&reference, &penalties, &cohort);
    n->state->initialize_probability_at_site(0, A);
    REQUIRE(correct_state.S == n->state->S);
    n = n->add_child_copying_state(C);
    correct_state.extend_probability_at_site(1, C);
    n->state->extend_probability_at_site(1, C);
    REQUIRE(correct_state.S == n->state->S);
  }
}


TEST_CASE( "Haplotype Manager can translate positions and identify and translate sites" , "[manager][indices] ") {
  string ref_refstring = "AAAAAAAAAAAAAAAAAAAA";
  string read_refstring = "AAAAAAAAAAAAAAAAAA";
  vector<size_t> ref_sites = {1,3,4,8,13,19};
  vector<size_t> read_sites = {3,5,6,13,17};
  linearReferenceStructure reference = build_ref(ref_refstring, ref_sites);
  vector<vector<alleleValue> > haplotypes = {
   {A, A, A, A, A, A}
  };
  haplotypeCohort cohort = haplotypeCohort(haplotypes, &reference);
  penaltySet penalties = penaltySet(-6, -9, haplotypes.size());
  referenceSequence ref_refseq = referenceSequence(ref_refstring);
  
  const char* ref_ref_cstring = ref_refstring.c_str();
  const char* read_ref_cstring = read_refstring.c_str();
  haplotypeManager hap_manager = haplotypeManager(
              &reference, 
              &cohort, 
              &penalties, 
              ref_ref_cstring,
              read_sites, 
              read_ref_cstring,
              2);
  SECTION( "Haplotype Manager can recognize its basic attributes" ) {
    REQUIRE(hap_manager.length() == 18);
    REQUIRE(hap_manager.read_sites() == 5);
    REQUIRE(hap_manager.shared_sites() == 2);
  }
  SECTION( "Haplotype manager can translate positions" ) {
    REQUIRE(hap_manager.ref_position(4) == 6);
    REQUIRE(hap_manager.read_position(3) == 1);
    
    REQUIRE(hap_manager.get_read_site_ref_position(0) == 5);
    REQUIRE(hap_manager.get_read_site_read_position(1) == 5);
    REQUIRE(hap_manager.get_read_site_ref_position(1) == 7);
    REQUIRE(hap_manager.read_index_to_shared_index(4) == 1);
    REQUIRE(hap_manager.read_index_to_read_only_index(3) == 2);
    REQUIRE(hap_manager.final_read_site_read_index() == 4);
    REQUIRE(hap_manager.final_read_site_read_position() == 17);
    REQUIRE(hap_manager.get_ref_site_below_read_site(2) == 3);
    REQUIRE(hap_manager.get_ref_site_below_read_site(3) == 4);
    
    REQUIRE(hap_manager.get_ref_site_ref_position(2) == 4);
    REQUIRE(hap_manager.final_ref_site() == 5);
    REQUIRE(hap_manager.final_ref_site_read_position() == 17);
    
    REQUIRE(hap_manager.shared_index_to_read_index(1) == 4);
    REQUIRE(hap_manager.shared_index_to_ref_index(1) == 5);
    REQUIRE(hap_manager.get_shared_site_ref_position(0) == 8);
    REQUIRE(hap_manager.get_shared_site_read_position(0) == 6);
    REQUIRE(hap_manager.get_shared_site_read_position(1) == 17);
  }
}

TEST_CASE( "Haplotype manager correctly handles tree-of-states", "[manager][tree]" ) {
  linearReferenceStructure reference = build_ref("AAAA", {0,1,2,3});
  vector<vector<alleleValue> > haplotypes = {
    {A, A, A, A},
    {A, T, A, A},
    {A, A, T, A},
    {A, A, A, T},
    {A, A, A, A}
  };
  haplotypeCohort cohort = haplotypeCohort(haplotypes, &reference);
  penaltySet penalties = penaltySet(-6, -9, haplotypes.size());
  string ref_string = "AAAA";
  string read_string = "AAAA";
  vector<size_t> read_sites = {0,1,2,3};
  
  SECTION( "HaplotypeManager can build a tree" ) {
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
    REQUIRE(true);
  }

  // SECTION( "HaplotypeManager maintains correct states for simple sets of sites" ) {
  //   haplotypeMatrix correct_state = haplotypeMatrix(
  //               &reference,
  //               &penalties,
  //               &cohort);
  //   
  //   haplotypeManager hap_manager = haplotypeManager(
  //             &reference,
  //             &cohort,
  //             &penalties,
  //             ref_string.c_str(),
  //             read_sites,
  //             read_string.c_str(),
  //             0);
  //     
  //   hap_manager.build_entire_tree(0);
  //   // hap_manager.print_tree();
  //   
  //   haplotypeStateNode* n = hap_manager.get_tree()->root;
  //   correct_state.initialize_probability_at_site(0, A);
  //   n = n->get_child(A);
  //   REQUIRE(correct_state.S == n->prefix_likelihood());
  //   correct_state.extend_probability_at_site(1, C);
  //   n = n->get_child(C);
  //   REQUIRE(correct_state.S == n->prefix_likelihood());
  //   correct_state.extend_probability_at_site(2, T);
  //   n = n->get_child(T);
  //   REQUIRE(correct_state.S == n->prefix_likelihood());
  //   correct_state.extend_probability_at_site(3, G);
  //   n = n->get_child(G);
  //   REQUIRE(correct_state.S == n->prefix_likelihood());
  // }
}

TEST_CASE( "Haplotype manager performs correct calculations in the presence of unshared sites", "[manager]" ) {
  SECTION( "All sites shared" ) {
    string ref_refstring = "AAAAAAA";
    string read_refstring = "AAAAAAA";
    vector<size_t> ref_sites = {0, 2, 5};
    vector<size_t> read_sites = {0, 2, 5};
    
    linearReferenceStructure reference = build_ref(ref_refstring, ref_sites);
    vector<vector<alleleValue> > haplotypes = {
      {A, A, A},
      {A, A, A},
      {A, T, A},
      {A, A, T},
      {A, G, C}
    };
    haplotypeCohort cohort = haplotypeCohort(haplotypes, &reference);
    penaltySet penalties = penaltySet(-6, -9, haplotypes.size());
    
    haplotypeMatrix AAA_state = haplotypeMatrix(&reference, &penalties, &cohort);
    haplotypeMatrix ATG_state = haplotypeMatrix(&reference, &penalties, &cohort);
    haplotypeMatrix TCG_state = haplotypeMatrix(&reference, &penalties, &cohort);
    vector<alleleValue> AAA = {A, A, A};
    vector<alleleValue> ATG = {A, T, G};
    vector<alleleValue> TCG = {T, C, G};
    
    AAA_state.initialize_probability_at_site(0, A);
    AAA_state.extend_probability_at_span_after((size_t)0, 0);
    AAA_state.extend_probability_at_site(1, A);
    AAA_state.extend_probability_at_span_after((size_t)1, 0);
    AAA_state.extend_probability_at_site(2, A);
    ATG_state.initialize_probability_at_site(0, A);
    ATG_state.extend_probability_at_span_after((size_t)0, 0);
    ATG_state.extend_probability_at_site(1, T);
    ATG_state.extend_probability_at_span_after((size_t)1, 0);
    ATG_state.extend_probability_at_site(2, G);
    TCG_state.initialize_probability_at_site(0, T);
    TCG_state.extend_probability_at_span_after((size_t)0, 0);
    TCG_state.extend_probability_at_site(1, C);
    TCG_state.extend_probability_at_span_after((size_t)1, 0);
    TCG_state.extend_probability_at_site(2, G);
    
    haplotypeManager hap_manager = haplotypeManager(
                &reference, 
                &cohort, 
                &penalties, 
                ref_refstring.c_str(),
                read_sites, 
                read_refstring.c_str(),
                0);
    hap_manager.build_entire_tree(0);
    hap_manager.print_tree();
    
    haplotypeStateNode* n;
    n = hap_manager.get_tree()->alleles_to_state(AAA);
    REQUIRE(n->prefix_likelihood() == AAA_state.S);
    n = hap_manager.get_tree()->alleles_to_state(ATG);
    REQUIRE(n->prefix_likelihood() == ATG_state.S);
    n = hap_manager.get_tree()->alleles_to_state(TCG);
    REQUIRE(n->prefix_likelihood() == TCG_state.S);
  }
  
  SECTION( "No sites" ) {
    string ref_refstring = "AAAAAAAAAAAA";
    string read_refstring = "AAAAAAA";
    vector<size_t> ref_sites = {0, 1, 9};
    vector<size_t> read_sites = {4};
    
    linearReferenceStructure reference = build_ref(ref_refstring, ref_sites);
    vector<vector<alleleValue> > haplotypes = {
      {A, A, A},
      {A, A, A},
      {A, T, A},
      {A, A, T},
      {A, G, C}
    };
    haplotypeCohort cohort = haplotypeCohort(haplotypes, &reference);
    penaltySet penalties = penaltySet(-6, -9, haplotypes.size());
    
    haplotypeManager hap_manager = haplotypeManager(
                &reference, 
                &cohort, 
                &penalties, 
                ref_refstring.c_str(),
                read_sites, 
                read_refstring.c_str(),
                2);
    hap_manager.build_entire_tree(0);
    // hap_manager.print_tree();
    
    REQUIRE(hap_manager.contains_ref_sites() == false);
    REQUIRE(hap_manager.length() == 7);
    REQUIRE(hap_manager.get_tree()->root->number_of_children() == 0);
  }
  
  SECTION( "Haplotype manager performs correct calculations without pruning") {
    string ref_refstring = "AAAAAAAAAAAAAAAAAAAA";
    string read_refstring = "AAAAAAAAAAAAAAAAAA";
    vector<size_t> ref_sites = {1,3,4,8,13,19};
    vector<size_t> read_sites = {3,5,6,13,17};
    linearReferenceStructure reference = build_ref(ref_refstring, ref_sites);
    vector<vector<alleleValue> > haplotypes = {
     {A, A, A, A, A, A},
     {C, T, A, A, A, A},
     {A, C, T, A, A, A},
     {A, A, C, T, A, A},
     {A, A, A, C, T, A},
     {A, A, A, A, C, T}
    };
    haplotypeCohort cohort = haplotypeCohort(haplotypes, &reference);
    penaltySet penalties = penaltySet(-6, -9, haplotypes.size());
    
    haplotypeMatrix twoAs_state = haplotypeMatrix(&reference, &penalties, &cohort);
    haplotypeMatrix twoTs_state = haplotypeMatrix(&reference, &penalties, &cohort);
    haplotypeMatrix twoGs_state = haplotypeMatrix(&reference, &penalties, &cohort);
    haplotypeMatrix AT_state = haplotypeMatrix(&reference, &penalties, &cohort);
    haplotypeMatrix root_state = haplotypeMatrix(&reference, &penalties, &cohort);
    
    twoAs_state.initialize_probability_at_span(1, 0);
    twoAs_state.extend_probability_at_site(1, A);
    twoAs_state.extend_probability_at_site(2, A);
    twoAs_state.extend_probability_at_span_after((size_t)2, 0);
    twoAs_state.extend_probability_at_site(3, A);
    twoAs_state.extend_probability_at_span_after((size_t)3, 0);
    twoAs_state.extend_probability_at_site(4, A);
    twoAs_state.extend_probability_at_span_after((size_t)4, 0);
    twoAs_state.extend_probability_at_site(5, A);
    
    twoTs_state.initialize_probability_at_span(1, 0);
    twoTs_state.extend_probability_at_site(1, A);
    twoTs_state.extend_probability_at_site(2, A);
    twoTs_state.extend_probability_at_span_after((size_t)2, 0);
    twoTs_state.extend_probability_at_site(3, T);
    twoTs_state.extend_probability_at_span_after((size_t)3, 0);
    twoTs_state.extend_probability_at_site(4, A);
    twoTs_state.extend_probability_at_span_after((size_t)4, 0);
    twoTs_state.extend_probability_at_site(5, T);
    
    twoGs_state.initialize_probability_at_span(1, 0);
    twoGs_state.extend_probability_at_site(1, A);
    twoGs_state.extend_probability_at_site(2, A);
    twoGs_state.extend_probability_at_span_after((size_t)2, 0);
    twoGs_state.extend_probability_at_site(3, G);
    twoGs_state.extend_probability_at_span_after((size_t)3, 0);
    twoGs_state.extend_probability_at_site(4, A);
    twoGs_state.extend_probability_at_span_after((size_t)4, 0);
    twoGs_state.extend_probability_at_site(5, G);
    
    AT_state.initialize_probability_at_span(1, 0);
    AT_state.extend_probability_at_site(1, A);
    AT_state.extend_probability_at_site(2, A);
    AT_state.extend_probability_at_span_after((size_t)2, 0);
    AT_state.extend_probability_at_site(3, A);
    AT_state.extend_probability_at_span_after((size_t)3, 0);
    AT_state.extend_probability_at_site(4, A);
    AT_state.extend_probability_at_span_after((size_t)4, 0);
    AT_state.extend_probability_at_site(5, T);
    
    root_state.initialize_probability_at_span(1, 0);
    root_state.extend_probability_at_site(1, A);
    root_state.extend_probability_at_site(2, A);
    root_state.extend_probability_at_span_after((size_t)2, 0);
  
    vector<alleleValue> twoAs = {A, A};
    vector<alleleValue> twoTs = {T, T};
    vector<alleleValue> twoGs = {G, G};
    vector<alleleValue> AT = {A, T};
  
    haplotypeStateNode* n;
    
    haplotypeManager hap_manager_root = haplotypeManager(
                &reference, 
                &cohort, 
                &penalties, 
                ref_refstring.c_str(),
                read_sites, 
                read_refstring.c_str(),
                2);
    haplotypeManager hap_manager = haplotypeManager(
                &reference, 
                &cohort, 
                &penalties, 
                ref_refstring.c_str(),
                read_sites, 
                read_refstring.c_str(),
                2);
    hap_manager_root.initialize_tree();
    n = hap_manager_root.get_tree()->root;
    // hap_manager.print_tree();
    REQUIRE(n->prefix_likelihood() == root_state.S);

    hap_manager.build_entire_tree(0);
    n = hap_manager.get_tree()->alleles_to_state(twoAs);
    REQUIRE(n->prefix_likelihood() == twoAs_state.S);
    n = hap_manager.get_tree()->alleles_to_state(twoTs);
    REQUIRE(n->prefix_likelihood() == twoTs_state.S);
    n = hap_manager.get_tree()->alleles_to_state(twoTs);
    REQUIRE(n->prefix_likelihood() == twoTs_state.S);
    n = hap_manager.get_tree()->alleles_to_state(AT);
    REQUIRE(n->prefix_likelihood() == AT_state.S);
  }
  
  /*
  SECTION( "Haplotype manager performs correct calculations with pruning") {
    haplotypeManager hap_manager = haplotypeManager(
                &reference, 
                &cohort, 
                &penalties, 
                ref_refstring.c_str(),
                read_sites, 
                read_refstring.c_str(),
                2);
    hap_manager.build_entire_tree(-24);
    haplotypeStateNode* n;
    n = hap_manager.get_tree()->alleles_to_state(sixAs);
    n = hap_manager.get_tree()->alleles_to_state(sixAs);
    REQUIRE(n->prefix_likelihood() == );
    n = hap_manager.get_tree()->alleles_to_state(sixTs);
    REQUIRE(n == nullptr);
    n = hap_manager.get_tree()->alleles_to_state(sixCs);
    REQUIRE(n == nullptr);
    n = hap_manager.get_tree()->alleles_to_state(ATTTTT);
    REQUIRE(n == nullptr);
    n = hap_manager.get_tree()->alleles_to_state(AAAATT);
    REQUIRE(n->prefix_likelihood() == );
  }*/
}
