#define CATCH_CONFIG_MAIN

#include <cmath>
#include "lh_reference.hpp"
#include "lh_probability.hpp"
#include "lh_input_haplotype.hpp"
#include "catch.hpp"

using namespace std;

linearReferenceStructure build_ref(string ref_seq, vector<size_t> positions) {
  vector<alleleValue> ref_values;
  for(size_t i = 0; i < positions.size(); i++) {
    ref_values.push_back(char_to_allele(ref_seq[positions[i]]));
  }
  return linearReferenceStructure(positions, ref_seq.length(), ref_values);
}

TEST_CASE( "linearReferenceStructure has desired structure", "[reference]" ) {
  // indices        0123456
  // sites           x  xx
  string ref_seq = "GATTACA";
  vector<size_t> positions = {1,4,5};
  linearReferenceStructure ref_struct = build_ref(ref_seq, positions);
  SECTION( "accessors give correct values" ) {
    REQUIRE(ref_struct.number_of_sites() == 3);
    REQUIRE(ref_struct.is_site(0) == false);
    REQUIRE(ref_struct.is_site(1) == true);
    REQUIRE(ref_struct.get_site_index(1) == 0);
    REQUIRE(ref_struct.get_site_index(5) == 2);
    REQUIRE(ref_struct.get_position(0) == 1);
    REQUIRE(ref_struct.get_position(2) == 5);
    REQUIRE(ref_struct.has_span_before(0) == true);
    REQUIRE(ref_struct.has_span_after(0) == true);
    REQUIRE(ref_struct.has_span_after(1) == false);
    REQUIRE(ref_struct.span_length_before(0) == 1);
    REQUIRE(ref_struct.span_length_after(0) == 2);
    REQUIRE(ref_struct.span_length_after(1) == 0);
    REQUIRE(ref_struct.span_length_after(2) == 1);
    REQUIRE(ref_struct.absolute_length() == 7);
  }
  SECTION( "build-from-strings gives same result as direct construction" ) {
    vector<string> haplotypes = {
      "GCTTA-A",
      "GATT-CA"
    };
    linearReferenceStructure ref_struct_from_strings =
              linearReferenceStructure(haplotypes, ref_seq);
    REQUIRE(ref_struct_from_strings.is_site(0) == false);
    REQUIRE(ref_struct_from_strings.is_site(1) == true);
    REQUIRE(ref_struct_from_strings.is_site(2) == false);
    REQUIRE(ref_struct_from_strings.is_site(3) == false);
    REQUIRE(ref_struct_from_strings.is_site(4) == true);
    REQUIRE(ref_struct_from_strings.is_site(5) == true);
    REQUIRE(ref_struct_from_strings.is_site(6) == false);
    REQUIRE(ref_struct_from_strings.get_site_index(1) == 0);
    REQUIRE(ref_struct_from_strings.get_site_index(5) == 2);
    REQUIRE(ref_struct_from_strings.get_position(0) == 1);
    REQUIRE(ref_struct_from_strings.get_position(2) == 5);
    REQUIRE(ref_struct_from_strings.span_length_before(0) == 1);
    REQUIRE(ref_struct_from_strings.span_length_after(0) == 2);
    REQUIRE(ref_struct_from_strings.span_length_after(1) == 0);
    REQUIRE(ref_struct_from_strings.span_length_after(2) == 1);
    REQUIRE(ref_struct_from_strings.has_span_before(0) == true);
    REQUIRE(ref_struct_from_strings.has_span_after(0) == true);
    REQUIRE(ref_struct_from_strings.has_span_after(1) == false);
  }
}

TEST_CASE( "haplotypeCohort construction behaves as desired", "[haplotype][reference]") {
  string ref_seq = "GATTACA";
  vector<size_t> positions = {1,4,5};
  linearReferenceStructure ref_struct = build_ref(ref_seq, positions);
  SECTION( "direct build works as intended" ) {
    // ref    GATTACA
    // sites   0  12
    // h_1    GATTAAA
    // h_2    GTTTAAA
    // h_3    GCTTAAA
    // h_4    GGTTAAA
    // h_5    G-TTAAA
    // h_6    G-TTAAA
    // As      1  66
    // Ts      1  00
    // Cs      1  00
    // Gs      1  00
    // -s      2  00
    // #       5  11
    vector<vector<alleleValue> > haplotypes = {
      {A,A,A},
      {T,A,A},
      {C,A,A},
      {G,A,A},
      {gap,A,A},
      {gap,A,A}
    };
    haplotypeCohort direct_cohort = haplotypeCohort(haplotypes, &ref_struct);
    
    REQUIRE(direct_cohort.size() == 6);
    REQUIRE(direct_cohort.allele_at(0,0) == A);
    REQUIRE(direct_cohort.allele_at(0,1) == T);
    REQUIRE(direct_cohort.allele_at(0,2) == C);
    REQUIRE(direct_cohort.allele_at(0,3) == G);
    REQUIRE(direct_cohort.allele_at(0,4) == gap);
    REQUIRE(direct_cohort.get_matches(0,A)[0] == 0);
    REQUIRE(direct_cohort.get_matches(0,T)[0] == 1);
    REQUIRE(direct_cohort.get_matches(0,C)[0] == 2);
    REQUIRE(direct_cohort.get_matches(0,G)[0] == 3);
    REQUIRE(direct_cohort.get_matches(0,gap)[0] == 4);
    REQUIRE(direct_cohort.get_matches(0,gap)[1] == 5);
    REQUIRE(direct_cohort.get_matches(1,A).size() == 6);
    REQUIRE(direct_cohort.number_matching(0,A) == 1);
    REQUIRE(direct_cohort.number_matching(0,T) == 1);
    REQUIRE(direct_cohort.number_matching(0,C) == 1);
    REQUIRE(direct_cohort.number_matching(0,G) == 1);
    REQUIRE(direct_cohort.number_matching(0,gap) == 2);
    REQUIRE(direct_cohort.number_matching(1,A) == 6);
    REQUIRE(direct_cohort.number_matching(1,T) == 0);
    REQUIRE(direct_cohort.number_not_matching(0,gap) == 4);
    REQUIRE(direct_cohort.number_not_matching(1,A) == 0);
  }
  SECTION( "string build works as intended" ) {
    // ref    GATTACA
    // sites   0  12
    // h_1    GATTAAA
    // h_2    GTTTAAA
    // h_3    GCTTAAA
    // h_4    GGTTAAA
    // h_5    G-TTAAA
    // h_6    G-TTAAA
    // As      1  66
    // Ts      1  00
    // Cs      1  00
    // Gs      1  00
    // -s      2  00
    // #       5  11
    vector<string> haplotypes {
      "GATTAAA",
      "GTTTAAA",
      "GCTTAAA",
      "GGTTAAA",
      "G-TTAAA",
      "G-TTAAA"
    };
    haplotypeCohort string_cohort = haplotypeCohort(haplotypes, &ref_struct);
    
    REQUIRE(string_cohort.size() == 6);
    REQUIRE(string_cohort.allele_at(0,0) == A);
    REQUIRE(string_cohort.allele_at(0,1) == T);
    REQUIRE(string_cohort.allele_at(0,2) == C);
    REQUIRE(string_cohort.allele_at(0,3) == G);
    REQUIRE(string_cohort.allele_at(0,4) == gap);
    REQUIRE(string_cohort.get_matches(0,A)[0] == 0);
    REQUIRE(string_cohort.get_matches(0,T)[0] == 1);
    REQUIRE(string_cohort.get_matches(0,C)[0] == 2);
    REQUIRE(string_cohort.get_matches(0,G)[0] == 3);
    REQUIRE(string_cohort.get_matches(0,gap)[0] == 4);
    REQUIRE(string_cohort.get_matches(0,gap)[1] == 5);
    REQUIRE(string_cohort.get_matches(1,A).size() == 6);
    REQUIRE(string_cohort.number_matching(0,A) == 1);
    REQUIRE(string_cohort.number_matching(0,T) == 1);
    REQUIRE(string_cohort.number_matching(0,C) == 1);
    REQUIRE(string_cohort.number_matching(0,G) == 1);
    REQUIRE(string_cohort.number_matching(0,gap) == 2);
    REQUIRE(string_cohort.number_matching(1,A) == 6);
    REQUIRE(string_cohort.number_matching(1,T) == 0);
    REQUIRE(string_cohort.number_not_matching(0,gap) == 4);
    REQUIRE(string_cohort.number_not_matching(1,A) == 0);
  }
}

TEST_CASE( "inputHaplotype methods behave as desired", "[haplotype]" ) {
  // indices        0123456
  // sites           x  xx
  string ref_seq = "GATTACA";
  vector<size_t> positions = {1,4,5};
  linearReferenceStructure ref_struct = build_ref(ref_seq, positions);
  
  SECTION( "augmentations are indexed as desired " ) {
    vector<alleleValue> alleles = {A, A};
    vector<size_t> counts = {1,2,3};
    inputHaplotype query = inputHaplotype(alleles, counts);
    REQUIRE(query.get_augmentations(-1) == 1);
    REQUIRE(query.get_augmentations(1) == 3);
  }
  SECTION( "build-from-string gives same result as direct construction" ) {
    // using ref_test from above
    vector<alleleValue> alleles = {A, A, C};
    // spans        o oo  o
    string seq_1 = "GATTACA";
    string seq_2 = "AAAAACC";
    vector<size_t> count_1 = {0,0,0,0};
    vector<size_t> count_2 = {1,2,0,1};
    inputHaplotype query_1_direct = inputHaplotype(alleles, count_1);
    inputHaplotype query_2_direct = inputHaplotype(alleles, count_2);
    inputHaplotype query_1_string = inputHaplotype(seq_1, ref_seq, &ref_struct);
    inputHaplotype query_2_string = inputHaplotype(seq_2, ref_seq, &ref_struct);
    
    bool initial_augs_match_zero = (query_1_direct.get_augmentations(-1) ==
              query_1_string.get_augmentations(-1));
    bool augs_match_zero = (query_1_direct.get_augmentations(0) == 
              query_1_string.get_augmentations(0));
    bool initial_augs_match_nonzero = (query_2_direct.get_augmentations(-1) == 
              query_2_string.get_augmentations(-1));
    bool augs_match_nonzero = (query_2_direct.get_augmentations(0) == 
                      query_2_string.get_augmentations(0));
    REQUIRE(initial_augs_match_zero);
    REQUIRE(augs_match_zero);
    REQUIRE(initial_augs_match_nonzero);
    REQUIRE(augs_match_nonzero);
    REQUIRE(query_1_direct.get_allele(0) == query_1_string.get_allele(0));
    REQUIRE(query_2_direct.get_allele(0) == query_2_string.get_allele(0));
  }
  SECTION( "edit-by-char gives same result as direct edit" ) {
    vector<alleleValue> alleles = {A, A, C};
    // indices      0123456
    // spans        o oo  o
    // ref_seq      GATTACA
    string seq_1 = "GATCACC";
    string seq_2 = "AATCACA";
    vector<size_t> count_1 = {0,1,0,1};
    vector<size_t> count_2 = {1,1,0,0};
    inputHaplotype query = inputHaplotype(alleles, count_1);
    inputHaplotype query_after = inputHaplotype(alleles, count_2);
    size_t i = 0;
    query.edit(i, seq_2[i], seq_1[i], ref_seq[i]);
    REQUIRE(query.get_augmentations(i) == query_after.get_augmentations(i));
    i = 2;
    query.edit(i, seq_2[i], seq_1[i], ref_seq[i]);
    REQUIRE(query.get_augmentations(i) == query_after.get_augmentations(i));
    i = 3;
    query.edit(i, seq_2[i], seq_1[i], ref_seq[i]);
    REQUIRE(query.get_augmentations(i) == query_after.get_augmentations(i));
    i = 6;
    query.edit(i, seq_2[i], seq_1[i], ref_seq[i]);
    REQUIRE(query.get_augmentations(i) == query_after.get_augmentations(i));
  }
}

TEST_CASE( "Haplotype probabilities are correctly calculated", "[haplotype][probability]" ) {
  // We need our penaltySet to work correctly for any likelihood calculations to
  double eps = 0.0000001;
  
  penaltySet pen_test = penaltySet(-6, -9, 10);
  // REQUIRE(fabs(pen_test.log1plusHminusRho - ) < eps);
  // REQUIRE(fabs(pen_test.log1minusRho - ) < eps);
  // REQUIRE(fabs(pen_test.log1minusMu - ) < eps);

  SECTION( "Partial likelihood is correctly calculated at an initial site" ) {
    string ref_seq = "A";
    vector<size_t> positions = {0};
    linearReferenceStructure ref_struct = build_ref(ref_seq, positions);
    
    vector<string> allA = vector<string>(2, "A");
    vector<string> allT = vector<string>(2, "T");
    vector<string> AT = {"A","T"};    
    
    haplotypeCohort cohort_allA = haplotypeCohort(allA, &ref_struct);
    haplotypeCohort cohort_allT = haplotypeCohort(allA, &ref_struct);
    haplotypeCohort cohort_AT = haplotypeCohort(allA, &ref_struct);
    
    inputHaplotype query = inputHaplotype((string)"A", ref_seq, &ref_struct);
    
    haplotypeMatrix matrix_allA = haplotypeMatrix(&ref_struct, &pen_test, 
                &cohort_allA, &query);
    double probability_allA = matrix_allA.calculate_probabilities();
    
    haplotypeMatrix matrix_allT = haplotypeMatrix(&ref_struct, &pen_test, 
                &cohort_allT, &query);
    double probability_allT = matrix_allT.calculate_probabilities();
    
    haplotypeMatrix matrix_AT = haplotypeMatrix(&ref_struct, &pen_test, 
                &cohort_AT, &query);
    double probability_AT = matrix_AT.calculate_probabilities();
    
    // TODO: correct values
    double R0_allA_p = 0;
    double R0_allT_p = 0;
    double R0_AT_p = 0;
    double R1_AT_p = 0;
    double p_allA_p = 0;
    double p_allT_p = 0;
    double p_AT_p = 0;
        
    bool R0_allA_correct = (fabs(matrix_allA.R[0][0] - R0_allA_p) < eps);
    
    bool R0_allT_correct = (fabs(matrix_allT.R[0][0] - R0_allT_p) < eps);
    
    bool R0_AT_correct = (fabs(matrix_AT.R[0][0] - R0_AT_p) < eps);
    bool R1_AT_correct = (fabs(matrix_AT.R[0][1] - R1_AT_p) < eps);
    
    bool probability_allA_correct = (fabs(probability_allA - p_allA_p) < eps);
    bool probability_allT_correct = (fabs(probability_allT - p_allT_p) < eps);
    bool probability_AT_correct = (fabs(probability_AT - p_AT_p) < eps);
    
    REQUIRE(R0_allA_correct);
    REQUIRE(R0_allT_correct);
    REQUIRE(R0_AT_correct);
    REQUIRE(R1_AT_correct);
    REQUIRE(probability_allA_correct);
    REQUIRE(probability_allT_correct);
    REQUIRE(probability_AT_correct);
  }
  SECTION( "Partial likelihood is correctly calculated at a noninitial site" ) {
    string ref_seq = "A";
    vector<size_t> positions = {0};
    linearReferenceStructure ref_struct = build_ref(ref_seq, positions);
    
    vector<string> allA = {"AA","TA","AA","TA"};
    vector<string> allT = {"AT","TT","AT","TT"};
    vector<string> AT = {"AA","TA","AT","TT"};    
    
    haplotypeCohort cohort_allA = haplotypeCohort(allA, &ref_struct);
    haplotypeCohort cohort_allT = haplotypeCohort(allA, &ref_struct);
    haplotypeCohort cohort_AT = haplotypeCohort(allA, &ref_struct);
    
    inputHaplotype query = inputHaplotype((string)"AA", ref_seq, &ref_struct);
    
    haplotypeMatrix matrix_allA = haplotypeMatrix(&ref_struct, &pen_test, 
                &cohort_allA, &query);
    double probability_allA = matrix_allA.calculate_probabilities();
    
    haplotypeMatrix matrix_allT = haplotypeMatrix(&ref_struct, &pen_test, 
                &cohort_allT, &query);
    double probability_allT = matrix_allT.calculate_probabilities();
    
    haplotypeMatrix matrix_AT = haplotypeMatrix(&ref_struct, &pen_test, 
                &cohort_AT, &query);
    double probability_AT = matrix_AT.calculate_probabilities();
    
    // TODO: fill in
    
    double R0_allA_p = 0;
    double R1_allA_p = 0;
    double R0_allT_p = 0;
    double R1_allT_p = 0;
    double R0_AT_p = 0;
    double R1_AT_p = 0;
    double R2_AT_p = 0;
    double R3_AT_p = 0;
    double p_allA_p = 0;
    double p_allT_p = 0;
    double p_AT_p = 0;    
    
    bool R0_allA_correct = (fabs(matrix_allA.R[1][0] - R0_allA_p) < eps);
    bool R1_allA_correct = (fabs(matrix_allA.R[1][1] - R0_allA_p) < eps);
    
    bool R0_allT_correct = (fabs(matrix_allT.R[1][0] - R0_allA_p) < eps);
    bool R1_allT_correct = (fabs(matrix_allT.R[1][0] - R0_allA_p) < eps);
    
    bool R0_AT_correct = (fabs(matrix_AT.R[1][0] - R0_AT_p) < eps);
    bool R1_AT_correct = (fabs(matrix_AT.R[1][1] - R0_AT_p) < eps);
    bool R2_AT_correct = (fabs(matrix_AT.R[1][2] - R0_AT_p) < eps);
    bool R3_AT_correct = (fabs(matrix_AT.R[1][3] - R0_AT_p) < eps);
    
    bool probability_allA_correct = (fabs(probability_allA - p_allA_p) < eps);
    bool probability_allT_correct = (fabs(probability_allT - p_allT_p) < eps);
    bool probability_AT_correct = (fabs(probability_AT - p_AT_p) < eps);
    
    REQUIRE(R0_allA_correct);
    REQUIRE(R1_allA_correct);
    REQUIRE(R0_allT_correct);
    REQUIRE(R1_allT_correct);
    REQUIRE(R0_AT_correct);
    REQUIRE(R1_AT_correct);
    REQUIRE(R2_AT_correct);
    REQUIRE(R3_AT_correct);
    REQUIRE(probability_allA_correct);
    REQUIRE(probability_allT_correct);
    REQUIRE(probability_AT_correct);
  }
  SECTION( "Partial likelihood is correctly calculated at a span following a site" ) {
    string ref_seq = "AAAAAA";
    vector<size_t> positions = {0};
    linearReferenceStructure ref_struct = build_ref(ref_seq, positions);
    
    vector<string> haplotypes = {
      "AAAAAA",
      "AAAAAA",
      "TAAAAA"
    };
    haplotypeCohort cohort = haplotypeCohort(haplotypes, &ref_struct);
    
    inputHaplotype query_0_aug = inputHaplotype((string)"AAAAAA", ref_seq, &ref_struct);
    inputHaplotype query_1_aug = inputHaplotype((string)"AAAAAT", ref_seq, &ref_struct);
    
    haplotypeMatrix matrix_0_aug = haplotypeMatrix(&ref_struct, &pen_test, 
                &cohort, &query_0_aug);
    double probability_0_aug = matrix_0_aug.calculate_probabilities();
    haplotypeMatrix matrix_1_aug = haplotypeMatrix(&ref_struct, &pen_test, 
                &cohort, &query_1_aug);
    double probability_1_aug = matrix_1_aug.calculate_probabilities();
    
    // TODO: fill in
    
    double p_0a_p = 0;
    double p_1a_p = 0;
    double R0_0a_p = 0;
    double R1_0a_p = 0;
    double R2_0a_p = 0;
    double R0_1a_p = 0;
    double R1_1a_p = 0;
    double R2_1a_p = 0;
    
    bool probability_0_aug_correct =
              (fabs(probability_0_aug - p_0a_p) < eps);
    bool probability_1_aug_correct =
              (fabs(probability_1_aug - p_0a_p) < eps);
    bool R0_0_aug_correct = (fabs(matrix_0_aug.R[0][0] - R0_0a_p) < eps);
    bool R1_0_aug_correct = (fabs(matrix_0_aug.R[0][1] - R1_0a_p) < eps);
    bool R2_0_aug_correct = (fabs(matrix_0_aug.R[0][2] - R2_0a_p) < eps);
    bool R0_1_aug_correct = (fabs(matrix_1_aug.R[0][0] - R0_1a_p) < eps);
    bool R1_1_aug_correct = (fabs(matrix_1_aug.R[0][1] - R1_1a_p) < eps);
    bool R2_1_aug_correct = (fabs(matrix_1_aug.R[0][2] - R2_1a_p) < eps);
    
    REQUIRE(R0_0_aug_correct);
    REQUIRE(R1_0_aug_correct);
    REQUIRE(R2_0_aug_correct);
    REQUIRE(R0_1_aug_correct);
    REQUIRE(R1_1_aug_correct);
    REQUIRE(R2_1_aug_correct);  
    REQUIRE(probability_1_aug_correct);
    REQUIRE(probability_1_aug_correct);
  }
  SECTION( "Partial likelihood is correctly calculated at an initial span" ) {
    string ref_seq = "AAAAAA";
    vector<size_t> positions = {5};
    linearReferenceStructure ref_struct = build_ref(ref_seq, positions);
    
    vector<string> haplotypes = {
      "AAAAAA",
      "AAAAAA",
      "AAAAAT"
    };
    haplotypeCohort cohort = haplotypeCohort(haplotypes, &ref_struct);
    
    inputHaplotype query_0_aug = inputHaplotype((string)"AAAAAA", ref_seq, &ref_struct);
    inputHaplotype query_1_aug = inputHaplotype((string)"TAAAAA", ref_seq, &ref_struct);
    
    haplotypeMatrix matrix_0_aug = haplotypeMatrix(&ref_struct, &pen_test, 
                &cohort, &query_0_aug);
    double probability_0_aug = matrix_0_aug.calculate_probabilities();
    haplotypeMatrix matrix_1_aug = haplotypeMatrix(&ref_struct, &pen_test, 
                &cohort, &query_1_aug);
    double probability_1_aug = matrix_1_aug.calculate_probabilities();
    
    // TODO: fill in
    
    double p_0a_p = 0;
    double p_1a_p = 0;
    double R0_0a_p = 0;
    double R1_0a_p = 0;
    double R2_0a_p = 0;
    double R0_1a_p = 0;
    double R1_1a_p = 0;
    double R2_1a_p = 0;
    
    bool probability_0_aug_correct =
              (fabs(probability_0_aug - p_0a_p) < eps);
    bool probability_1_aug_correct =
              (fabs(probability_1_aug - p_0a_p) < eps);
    bool R0_0_aug_correct = (fabs(matrix_0_aug.R[0][0] - R0_0a_p) < eps);
    bool R1_0_aug_correct = (fabs(matrix_0_aug.R[0][1] - R1_0a_p) < eps);
    bool R2_0_aug_correct = (fabs(matrix_0_aug.R[0][2] - R2_0a_p) < eps);
    bool R0_1_aug_correct = (fabs(matrix_1_aug.R[0][0] - R0_1a_p) < eps);
    bool R1_1_aug_correct = (fabs(matrix_1_aug.R[0][1] - R1_1a_p) < eps);
    bool R2_1_aug_correct = (fabs(matrix_1_aug.R[0][2] - R2_1a_p) < eps);
    
    REQUIRE(R0_0_aug_correct);
    REQUIRE(R1_0_aug_correct);
    REQUIRE(R2_0_aug_correct);
    REQUIRE(R0_1_aug_correct);
    REQUIRE(R1_1_aug_correct);
    REQUIRE(R2_1_aug_correct);  
    REQUIRE(probability_1_aug_correct);
    REQUIRE(probability_1_aug_correct);
  }
  SECTION( "Partial likelihood is correctly calculated at a series of sites" ) {
    string ref_seq = "AAAAA";
    vector<size_t> positions = {0,1,2,3,4};
    linearReferenceStructure ref_struct = build_ref(ref_seq, positions);
    
    vector<string> haplotypes = {
      "AAAAA",
      "TATAT",
      "ATACA"
    };
    // 10110
    // 00011
    // 11100
    haplotypeCohort cohort = haplotypeCohort(haplotypes, &ref_struct);
    
    string query_string = "ATAAT";
    inputHaplotype query = inputHaplotype(query_string, ref_seq, &ref_struct);
    
    haplotypeMatrix matrix = haplotypeMatrix(&ref_struct, &pen_test, 
                &cohort, &query);
    double probability = matrix.calculate_probabilities();
    
    //TODO: fill in
    
    vector<vector<double> > correct_Rs = {
      {0, 0, 0},
      {0, 0, 0},
      {0, 0, 0},
      {0, 0, 0},
      {0, 0, 0}
    };
    
    vector<double> correctSs = {0, 0, 0, 0, 0};
    
    bool R00_same = (fabs(matrix.R[0][0] == correct_Rs[0][0]) < eps);
    bool R01_same = (fabs(matrix.R[0][1] == correct_Rs[0][1]) < eps); 
    bool R02_same = (fabs(matrix.R[0][2] == correct_Rs[0][2]) < eps); 
    
    bool R10_same = (fabs(matrix.R[1][0] == correct_Rs[1][0]) < eps);
    bool R11_same = (fabs(matrix.R[1][1] == correct_Rs[1][1]) < eps); 
    bool R12_same = (fabs(matrix.R[1][2] == correct_Rs[1][2]) < eps); 
    
    bool R20_same = (fabs(matrix.R[2][0] == correct_Rs[2][0]) < eps);
    bool R21_same = (fabs(matrix.R[2][1] == correct_Rs[2][1]) < eps); 
    bool R22_same = (fabs(matrix.R[2][2] == correct_Rs[2][2]) < eps); 
    
    bool R30_same = (fabs(matrix.R[3][0] == correct_Rs[3][0]) < eps);
    bool R31_same = (fabs(matrix.R[3][1] == correct_Rs[3][1]) < eps); 
    bool R32_same = (fabs(matrix.R[3][2] == correct_Rs[3][2]) < eps); 
    
    bool R40_same = (fabs(matrix.R[4][0] == correct_Rs[4][0]) < eps);
    bool R41_same = (fabs(matrix.R[4][1] == correct_Rs[4][1]) < eps); 
    bool R42_same = (fabs(matrix.R[4][2] == correct_Rs[4][2]) < eps);
    
    bool S0_same = (fabs(matrix.S[0] == correctSs[0]) < eps);
    bool S1_same = (fabs(matrix.S[1] == correctSs[1]) < eps);
    bool S2_same = (fabs(matrix.S[2] == correctSs[2]) < eps);
    bool S3_same = (fabs(matrix.S[3] == correctSs[3]) < eps);
    bool probabilities_same = (fabs(probability == correctSs[4]) < eps);
    
    REQUIRE(R00_same);
    REQUIRE(R01_same);
    REQUIRE(R02_same);
    REQUIRE(R10_same);
    REQUIRE(R11_same);
    REQUIRE(R12_same);
    REQUIRE(R20_same);
    REQUIRE(R21_same);
    REQUIRE(R22_same);
    REQUIRE(R30_same);
    REQUIRE(R31_same);
    REQUIRE(R32_same);
    REQUIRE(R40_same);
    REQUIRE(R41_same);
    REQUIRE(R42_same);
    REQUIRE(S0_same);
    REQUIRE(S1_same);
    REQUIRE(S2_same);
    REQUIRE(S3_same);
    REQUIRE(probabilities_same);
  }
  
  SECTION( "Partial likelihood at a series of sites without variation equals the likelihood at a span of equivalent length" ) {
    string ref_seq = "AAAAA";
    vector<size_t> positions_span = {0};
    linearReferenceStructure ref_struct_span = build_ref(ref_seq, 
                positions_span);

    vector<size_t> positions = {0,1,2,3,4};
    linearReferenceStructure ref_struct = build_ref(ref_seq, positions);
        
    vector<string> haplotypes = {
      "AAAAA",
      "AAAAA",
      "TAAAA"
    };
    haplotypeCohort cohort = haplotypeCohort(haplotypes, &ref_struct);
    haplotypeCohort cohort_span = haplotypeCohort(haplotypes, &ref_struct_span);    
    
    inputHaplotype query = inputHaplotype(ref_seq, ref_seq, &ref_struct);
    inputHaplotype query_span = inputHaplotype(ref_seq, ref_seq, 
                &ref_struct_span);
    
    haplotypeMatrix matrix = haplotypeMatrix(&ref_struct, &pen_test, 
                &cohort, &query);
    haplotypeMatrix matrix_span = haplotypeMatrix(&ref_struct_span, &pen_test, 
                &cohort_span, &query_span);
    double probability = matrix.calculate_probabilities();
    double probability_span = matrix_span.calculate_probabilities();
    
    bool S0_same = (fabs(matrix.S[0] - matrix_span.S[0]) < eps);
    bool probabilities_same = (fabs(probability - probability_span) < eps);
    
    bool R00_same = (fabs(matrix.R[0][0] - matrix_span.R[0][0]) < eps);
    bool R01_same = (fabs(matrix.R[0][1] - matrix_span.R[0][1]) < eps); 
    bool R02_same = (fabs(matrix.R[0][2] - matrix_span.R[0][2]) < eps); 
    
    bool R10_same = (fabs(matrix.R[4][1] - matrix_span.R[1][1]) < eps);
    bool R11_same = (fabs(matrix.R[4][1] - matrix_span.R[1][1]) < eps); 
    bool R12_same = (fabs(matrix.R[4][1] - matrix_span.R[1][1]) < eps); 
    
    REQUIRE(R00_same);
    REQUIRE(R01_same);
    REQUIRE(R02_same);
    REQUIRE(S0_same);
    REQUIRE(R10_same);
    REQUIRE(R11_same);
    REQUIRE(R12_same);
    REQUIRE(probabilities_same);
  }  
}