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
    inputHaplotype query_0_direct = inputHaplotype(alleles, count_1);
    inputHaplotype query_1_direct = inputHaplotype(alleles, count_2);
    inputHaplotype query_0_string = inputHaplotype(seq_1, ref_seq, &ref_struct);
    inputHaplotype query_1_string = inputHaplotype(seq_2, ref_seq, &ref_struct);
    
    REQUIRE(query_0_direct.get_augmentations(-1) ==
              query_0_string.get_augmentations(-1));
    REQUIRE(query_0_direct.get_augmentations(0) == 
              query_0_string.get_augmentations(0));
    REQUIRE(query_1_direct.get_augmentations(-1) == 
              query_1_string.get_augmentations(-1));
    REQUIRE(query_1_direct.get_augmentations(0) == 
              query_1_string.get_augmentations(0));
    REQUIRE(query_0_direct.get_allele(0) == query_0_string.get_allele(0));
    REQUIRE(query_1_direct.get_allele(0) == query_1_string.get_allele(0));
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
    inputHaplotype query = inputHaplotype(alleles, count_1, &ref_struct);
    inputHaplotype query_after = inputHaplotype(alleles, count_2, &ref_struct);
    size_t i = 0;
    query.edit(i, seq_2[i], seq_1[i], ref_seq[i]);
    REQUIRE(query.get_augmentations(0) == query_after.get_augmentations(0));
    i = 2;
    query.edit(i, seq_2[i], seq_1[i], ref_seq[i]);
    REQUIRE(query.get_augmentations(1) == query_after.get_augmentations(1));
    i = 3;
    query.edit(i, seq_2[i], seq_1[i], ref_seq[i]);
    REQUIRE(query.get_augmentations(1) == query_after.get_augmentations(1));
    i = 6;
    query.edit(i, seq_2[i], seq_1[i], ref_seq[i]);
    REQUIRE(query.get_augmentations(2) == query_after.get_augmentations(2));
  }
}

TEST_CASE( "Haplotype probabilities are correctly calculated", "[haplotype][probability]" ) {
  // We need our penaltySet to work correctly for any likelihood calculations to
  double eps = 0.0000001;
  
  SECTION( "Partial likelihood is correctly calculated at an initial site" ) {
    penaltySet penalties = penaltySet(-6, -9, 2);
    
    string ref_seq = "A";
    vector<size_t> positions = {0};
    linearReferenceStructure ref_struct = build_ref(ref_seq, positions);

    vector<string> allA = vector<string>(2, "A");
    vector<string> allT = vector<string>(2, "T");
    vector<string> AT = {"A","T"};    
    
    haplotypeCohort cohort_allA = haplotypeCohort(allA, &ref_struct);
    haplotypeCohort cohort_allT = haplotypeCohort(allT, &ref_struct);
    haplotypeCohort cohort_AT = haplotypeCohort(AT, &ref_struct);
    
    inputHaplotype query = inputHaplotype((string)"A", ref_seq, &ref_struct);

    haplotypeMatrix matrix_allA = haplotypeMatrix(&ref_struct, &penalties, 
                &cohort_allA, &query);
    double probability_allA = matrix_allA.calculate_probability();

    haplotypeMatrix matrix_allT = haplotypeMatrix(&ref_struct, &penalties, 
                &cohort_allT, &query);
    double probability_allT = matrix_allT.calculate_probability();

    haplotypeMatrix matrix_AT = haplotypeMatrix(&ref_struct, &penalties, 
                &cohort_AT, &query);
    double probability_AT = matrix_AT.calculate_probability();

    // TODO: correct values
    double R0_allA_p = log(0.5) + penalties.log_mu_complement;
    double R0_allT_p = log(0.5) + penalties.log_mu;
    double R0_AT_p = R0_allA_p;
    double R1_AT_p = R0_allT_p;
    double p_allA_p = log(2) + R0_allA_p;
    double p_allT_p = log(2) + R0_allT_p;
    double p_AT_p = logsum(R0_AT_p,R1_AT_p);
        
    REQUIRE(fabs(matrix_allA.R[0][0] - R0_allA_p) < eps);
    REQUIRE(fabs(matrix_allT.R[0][0] - R0_allT_p) < eps);
    
    REQUIRE(fabs(matrix_AT.R[0][0] - R0_AT_p) < eps);
    REQUIRE(fabs(matrix_AT.R[0][1] - R1_AT_p) < eps);
    
    REQUIRE(fabs(probability_allA - p_allA_p) < eps);
    REQUIRE(fabs(probability_allT - p_allT_p) < eps);
    REQUIRE(fabs(probability_AT - p_AT_p) < eps);
  }
  SECTION( "Partial likelihood is correctly calculated at a noninitial site" ) {
    penaltySet penalties = penaltySet(-6, -9, 4);
    
    string ref_seq = "AA";
    vector<size_t> positions = {0,1};
    linearReferenceStructure ref_struct = build_ref(ref_seq, positions);
    
    vector<string> allA = {"AA","TA","AA","TA"};
    vector<string> allT = {"AT","TT","AT","TT"};
    vector<string> AT = {"AA","TA","AT","TT"};    
    
    haplotypeCohort cohort_allA = haplotypeCohort(allA, &ref_struct);
    haplotypeCohort cohort_allT = haplotypeCohort(allT, &ref_struct);
    haplotypeCohort cohort_AT = haplotypeCohort(AT, &ref_struct);
    
    inputHaplotype query = inputHaplotype((string)"AA", ref_seq, &ref_struct);
    
    haplotypeMatrix matrix_allA = haplotypeMatrix(&ref_struct, &penalties, 
                &cohort_allA, &query);
    double probability_allA = matrix_allA.calculate_probability();
    
    haplotypeMatrix matrix_allT = haplotypeMatrix(&ref_struct, &penalties, 
                &cohort_allT, &query);
    double probability_allT = matrix_allT.calculate_probability();
    
    haplotypeMatrix matrix_AT = haplotypeMatrix(&ref_struct, &penalties, 
                &cohort_AT, &query);
    double probability_AT = matrix_AT.calculate_probability();
        
    double R0A = log(0.25) + penalties.log_mu_complement;
    double R0T = log(0.25) + penalties.log_mu;
    double S0 = log(0.5);
    double ftR0A = penalties.log_ft_base + R0A;
    double ftR0T = penalties.log_ft_base + R0T;
    double pS = penalties.log_rho + S0;
    
    double R0_allA_p = penalties.log_mu_complement + logsum(ftR0A, pS);
    double R1_allA_p = penalties.log_mu_complement + logsum(ftR0T, pS);
    double R0_allT_p = penalties.log_mu + logsum(ftR0A, pS);
    double R1_allT_p = penalties.log_mu + logsum(ftR0T, pS);
    double R0_AT_p = R0_allA_p;
    double R1_AT_p = R1_allA_p;
    double R2_AT_p = R0_allT_p;
    double R3_AT_p = R1_allT_p;
    double p_allA_p = log(0.5) + penalties.log_mu_complement + log1p(2*exp(-6));
    double p_allT_p = log(0.5) + penalties.log_mu + log1p(2*exp(-6));
    double p_AT_p = logsum(log(0.25) + penalties.log_ft_base, penalties.log_rho);    
    
    REQUIRE(fabs(matrix_allA.R[1][0] - R0_allA_p) < eps);
    REQUIRE(fabs(matrix_allA.R[1][1] - R1_allA_p) < eps);
    
    REQUIRE(fabs(matrix_allT.R[1][0] - R0_allT_p) < eps);
    REQUIRE(fabs(matrix_allT.R[1][1] - R1_allT_p) < eps);
    
    REQUIRE(fabs(matrix_AT.R[1][0] - R0_AT_p) < eps);
    REQUIRE(fabs(matrix_AT.R[1][1] - R1_AT_p) < eps);
    REQUIRE(fabs(matrix_AT.R[1][2] - R2_AT_p) < eps);
    REQUIRE(fabs(matrix_AT.R[1][3] - R3_AT_p) < eps);
    
    REQUIRE(fabs(probability_allA - p_allA_p) < eps);
    REQUIRE(fabs(probability_allT - p_allT_p) < eps);
    REQUIRE(fabs(probability_AT - p_AT_p) < eps);
  }
  SECTION( "Partial likelihood is correctly calculated at a span following a site" ) {
    penaltySet penalties = penaltySet(-6, -9, 3);
    
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
    
    haplotypeMatrix matrix_0_aug = haplotypeMatrix(&ref_struct, &penalties, 
                &cohort, &query_0_aug);
    double probability_0_aug = matrix_0_aug.calculate_probability();
    haplotypeMatrix matrix_1_aug = haplotypeMatrix(&ref_struct, &penalties, 
                &cohort, &query_1_aug);
    double probability_1_aug = matrix_1_aug.calculate_probability();
    
    double lmu = penalties.log_mu;
    double lmu_c = penalties.log_mu_complement;
    double lfsb = penalties.log_fs_base;
    double lftb = penalties.log_ft_base;
    double lS_i = logdiff(log(2), lmu) - log(3);
    double lrho = penalties.log_rho;
    
    double RHS = logsum(lrho + lS_i + 4*lfsb, 
                        lftb + lS_i + logdiff(4*lfsb, 4*lftb) - log(3));
    double lR_i = lmu_c - log(3);
    double lR_im = lmu - log(3);
    
    double R0_0a_p = lmu_c*5 + logsum(lftb*5 + lR_i, RHS);
    double R2_0a_p = lmu_c*5 + logsum(lftb*5 + lR_im, RHS);
    double R0_1a_p = lmu_c*4 + lmu + logsum(lftb*5 + lR_i, RHS);
    double R2_1a_p = lmu_c*4 + lmu + logsum(lftb*5 + lR_im, RHS);
    double p_0a_p = lmu_c*5 + logsum(lftb*5 + lS_i, RHS + log(3));
    double p_1a_p = lmu_c*4 + lmu + logsum(lftb*5 + lS_i, RHS + log(3));
    
    REQUIRE(fabs(matrix_0_aug.R[0][0] - R0_0a_p) < eps);
    REQUIRE(fabs(matrix_0_aug.R[0][2] - R2_0a_p) < eps);
    REQUIRE(fabs(matrix_1_aug.R[0][0] - R0_1a_p) < eps);
    REQUIRE(fabs(matrix_1_aug.R[0][2] - R2_1a_p) < eps);  
    REQUIRE(fabs(probability_0_aug - p_0a_p) < eps);
    REQUIRE(fabs(probability_1_aug - p_1a_p) < eps);
  }
  SECTION( "Partial likelihood is correctly calculated at an initial span" ) {
    penaltySet penalties = penaltySet(-6, -9, 3);
    
    string ref_seq = "AAAAAA";
    vector<size_t> positions = {5};
    linearReferenceStructure ref_struct = build_ref(ref_seq, positions);
    
    vector<string> haplotypes = {
      "AAAAAA",
      "AAAAAA",
      "AAAAAT"
    };
    haplotypeCohort cohort = haplotypeCohort(haplotypes, &ref_struct);
    
    inputHaplotype query_0_aug = 
              inputHaplotype((string)"AAAAAA", ref_seq, &ref_struct);
    inputHaplotype query_1_aug = 
              inputHaplotype((string)"TAAAAA", ref_seq, &ref_struct);
    
    haplotypeMatrix matrix_0_aug = haplotypeMatrix(&ref_struct, &penalties, 
                &cohort, &query_0_aug);
    double probability_0_aug = matrix_0_aug.calculate_probability();
    haplotypeMatrix matrix_1_aug = haplotypeMatrix(&ref_struct, &penalties, 
                &cohort, &query_1_aug);
    double probability_1_aug = matrix_1_aug.calculate_probability();
    
    double lmu = penalties.log_mu;
    double lmu_c = penalties.log_mu_complement;
    double lfsb = penalties.log_fs_base;
    
    double R0_0a_p = 6*lmu_c + 5*lfsb - log(3);
    double R2_0a_p = 5*lmu_c + 5*lfsb + lmu - log(3);
    double R0_1a_p = 5*lmu_c + lmu + 5*lfsb - log(3);
    double R2_1a_p = 4*lmu_c + 2*lmu + 5*lfsb - log(3);
    double p_0a_p = 5*lmu_c + 5*lfsb + logdiff(log(2), lmu) - log(3);
    double p_1a_p = 4*lmu_c + lmu + 5*lfsb + logdiff(log(2), lmu) - log(3);
    
    REQUIRE(fabs(matrix_0_aug.R[0][0] - R0_0a_p) < eps);
    REQUIRE(fabs(matrix_0_aug.R[0][2] - R2_0a_p) < eps);
    REQUIRE(fabs(matrix_1_aug.R[0][0] - R0_1a_p) < eps);
    REQUIRE(fabs(matrix_1_aug.R[0][2] - R2_1a_p) < eps);
    REQUIRE(fabs(probability_0_aug - p_0a_p) < eps);
    REQUIRE(fabs(probability_1_aug - p_1a_p) < eps);
  }
  SECTION( "Partial likelihood is correctly calculated at a series of sites" ) {
    penaltySet penalties = penaltySet(-6, -9, 3);
    
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
    
    haplotypeMatrix matrix = haplotypeMatrix(&ref_struct, &penalties, 
                &cohort, &query);
    double probability = matrix.calculate_probability();
    
    double lmu_c = penalties.log_mu_complement;
    double lmu = penalties.log_mu;
    double lrho = penalties.log_rho;
    double lftb = penalties.log_ft_base;
    double l2_mu = logdiff(log(2), lmu);
    
    vector<double> correct_Ss;
    vector<vector<double> > correct_Rs;
    
    // R0
    correct_Rs.push_back({
        -log(3) + lmu_c,
        -log(3) + lmu,
        -log(3) + lmu_c
    });
    // S0
    correct_Ss.push_back(-log(3) + l2_mu);
    // R1
    correct_Rs.push_back({
        lmu + logsum(lftb + correct_Rs[0][0], lrho + correct_Ss[0]),
        lmu + logsum(lftb + correct_Rs[0][1], lrho + correct_Ss[0]),
        lmu_c + logsum(lftb + correct_Rs[0][2], lrho + correct_Ss[0])
    });
    // S1
    correct_Ss.push_back(logsum(correct_Rs[1][0],
                                logsum(correct_Rs[1][1],correct_Rs[1][2])));
    // R2
    correct_Rs.push_back({
        lmu_c + logsum(lftb + correct_Rs[1][0], lrho + correct_Ss[1]),
        lmu + logsum(lftb + correct_Rs[1][1], lrho + correct_Ss[1]),
        lmu_c + logsum(lftb + correct_Rs[1][2], lrho + correct_Ss[1])
    });
    // S2
    correct_Ss.push_back(logsum(correct_Rs[2][0],
                                logsum(correct_Rs[2][1],correct_Rs[2][2])));
    // R3
    correct_Rs.push_back({
        lmu_c + logsum(lftb + correct_Rs[2][0], lrho + correct_Ss[2]),
        lmu_c + logsum(lftb + correct_Rs[2][1], lrho + correct_Ss[2]),
        lmu + logsum(lftb + correct_Rs[2][2], lrho + correct_Ss[2])
    });
    // S3
    correct_Ss.push_back(logsum(correct_Rs[3][0],
                                logsum(correct_Rs[3][1],correct_Rs[3][2])));
    // R4
    correct_Rs.push_back({
        lmu + logsum(lftb + correct_Rs[3][0], lrho + correct_Ss[3]),
        lmu_c + logsum(lftb + correct_Rs[3][1], lrho + correct_Ss[3]),
        lmu + logsum(lftb + correct_Rs[3][2], lrho + correct_Ss[3])
    });
    // S4
    correct_Ss.push_back(logsum(correct_Rs[4][0],
                                logsum(correct_Rs[4][1],correct_Rs[4][2])));
    
    REQUIRE(fabs(matrix.R[0][0] - correct_Rs[0][0]) < eps);
    REQUIRE(fabs(matrix.R[0][1] - correct_Rs[0][1]) < eps); 
    REQUIRE(fabs(matrix.R[0][2] - correct_Rs[0][2]) < eps); 
    REQUIRE(fabs(matrix.S[0] - correct_Ss[0]) < eps);

    REQUIRE(fabs(matrix.R[1][0] - correct_Rs[1][0]) < eps);
    REQUIRE(fabs(matrix.R[1][1] - correct_Rs[1][1]) < eps); 
    REQUIRE(fabs(matrix.R[1][2] - correct_Rs[1][2]) < eps); 
    REQUIRE(fabs(matrix.S[1] - correct_Ss[1]) < eps);

    REQUIRE(fabs(matrix.R[2][0] - correct_Rs[2][0]) < eps);
    REQUIRE(fabs(matrix.R[2][1] - correct_Rs[2][1]) < eps); 
    REQUIRE(fabs(matrix.R[2][2] - correct_Rs[2][2]) < eps); 
    REQUIRE(fabs(matrix.S[2] - correct_Ss[2]) < eps);

    REQUIRE(fabs(matrix.R[3][0] - correct_Rs[3][0]) < eps);
    REQUIRE(fabs(matrix.R[3][1] - correct_Rs[3][1]) < eps); 
    REQUIRE(fabs(matrix.R[3][2] - correct_Rs[3][2]) < eps); 
    REQUIRE(fabs(matrix.S[3] - correct_Ss[3]) < eps);

    REQUIRE(fabs(matrix.R[4][0] - correct_Rs[4][0]) < eps);
    REQUIRE(fabs(matrix.R[4][1] - correct_Rs[4][1]) < eps); 
    REQUIRE(fabs(matrix.R[4][2] - correct_Rs[4][2]) < eps);
    
    REQUIRE(fabs(probability - correct_Ss[4]) < eps);
  }
  
  SECTION( "Partial likelihood at a series of sites without variation equals the likelihood at a span of equivalent length" ) {
    penaltySet penalties = penaltySet(-6, -9, 3);
    
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
                
    haplotypeMatrix matrix = haplotypeMatrix(&ref_struct, &penalties, 
                &cohort, &query);
    haplotypeMatrix matrix_span = haplotypeMatrix(&ref_struct_span, &penalties, 
                &cohort_span, &query_span);
                
    double probability = matrix.calculate_probability();
    double probability_span = matrix_span.calculate_probability();
      
    REQUIRE(fabs(matrix.R[4][0] - matrix_span.R[0][0]) < eps);
    REQUIRE(fabs(matrix.R[4][1] - matrix_span.R[0][1]) < eps); 
    REQUIRE(fabs(matrix.R[4][2] - matrix_span.R[0][2]) < eps); 
    REQUIRE(fabs(probability - probability_span) < eps);
  }  
}

TEST_CASE( "Relative indexing works", "[haplotype][reference][input]" ) {
  //                01234567890123456789
  // sites              4    9    4
  string ref_seq = "AAAAAAAAAAAAAAAAAAAA";
  // left            AAA
  // middle              AAA 
  // right                         AAA    
  // 9 14                   AAAAAAAA
  //                01234567890123456789
  
  vector<size_t> positions = {4, 9, 14};
  linearReferenceStructure ref_struct = build_ref(ref_seq, positions);
  
  SECTION( "relative indices are produced as desired " ) {
    string eight = "TTTTTTTT";
    string three = "AAA";
    string six = "TTTTTT";
    inputHaplotype nosites_left = inputHaplotype(three, ref_seq, 
                                                 &ref_struct, 1, 3);
    inputHaplotype nosites_middle = inputHaplotype(three, ref_seq,
                                                  &ref_struct, 5, 3);
    inputHaplotype nosites_right = inputHaplotype(three, ref_seq,
                                                    &ref_struct, 15, 3);
    inputHaplotype overlap_9_16 = inputHaplotype(eight, ref_seq,
                                                 &ref_struct, 8, 8);
    inputHaplotype overlap_4_9 = inputHaplotype(six, ref_seq,
                                                 &ref_struct, 4, 6);
    REQUIRE(nosites_left.has_sites() == false);
    REQUIRE(nosites_right.has_sites() == false);
    REQUIRE(nosites_middle.has_sites() == false);
    REQUIRE(overlap_9_16.has_sites() == true);
    REQUIRE(overlap_9_16.number_of_sites() == 2);
    REQUIRE(overlap_9_16.get_left_tail() == 1);
    REQUIRE(overlap_9_16.get_rel_index(0) == 1);
    REQUIRE(overlap_9_16.get_augmentations(-1) == 1);
    REQUIRE(overlap_9_16.get_augmentations(0) == 4);
    REQUIRE(overlap_9_16.get_augmentations(1) == 1);
    REQUIRE(overlap_4_9.has_left_tail() == false);
    REQUIRE(overlap_4_9.number_of_sites() == 2);
    REQUIRE(overlap_4_9.has_span_after(1) == false);
    REQUIRE(overlap_9_16.get_augmentations(0) == 4);
  }
}