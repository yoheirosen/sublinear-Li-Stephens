#define CATCH_CONFIG_MAIN

#include <cmath>
#include "math.hpp"
#include "reference.hpp"
#include "probability.hpp"
#include "input_haplotype.hpp"
#include "delay_multiplier.hpp"
#include "catch.hpp"
#include <iostream>
#include <fstream>
#include <cstdio>

using namespace std;

siteIndex build_ref(const string& ref_seq, vector<size_t>& positions) {
  return siteIndex(positions, ref_seq.length());
}

TEST_CASE( "siteIndex structure and accessors", "[siteIndex]" ) {
  // indices        0123456
  // sites           x  xx
  string ref_seq = "GATTACA";
  vector<size_t> positions = {1,4,5};
  siteIndex ref_struct = build_ref(ref_seq, positions);
  SECTION( "siteIndex accessor correctness" ) {
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
    REQUIRE(ref_struct.length_in_bp() == 7);
  }
  SECTION( "build-from-strings finds all sites" ) {
    vector<string> haplotypes = {
      "GCTTA-A",
      "GATT-CA"
    };
    siteIndex ref_struct_from_strings = siteIndex(haplotypes);
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
  SECTION( "reading and writing to file works" ) {
    ofstream testout;
    testout.open("testout.linref", ios::out | ios::trunc);
    ref_struct.serialize_human(testout);
    testout.close();
    ifstream testin;
    testin.open("testout.linref", ios::in);
    siteIndex read_ref_struct(testin);
    testin.close();
    remove("testout.linref");
    REQUIRE(read_ref_struct.number_of_sites() == 3);
    REQUIRE(read_ref_struct.get_site_index(5) == 2);
    REQUIRE(read_ref_struct.length_in_bp() == 7);
    REQUIRE(read_ref_struct.span_length_after(2) == 1);
  }
}

TEST_CASE( "haplotypeCohort accessors", "[cohort][cohort-accessors]") {
  string ref_seq = "GATTACA";
  vector<size_t> positions = {1,4,5};
  siteIndex ref_struct = build_ref(ref_seq, positions);
  vector<vector<alleleValue> > haplotypes = {
    {A,A,A},
    {T,A,A},
    {C,A,T},
    {G,A,G},
    {gap,A,A},
    {gap,A,A}
  };
  haplotypeCohort cohort = haplotypeCohort(haplotypes, &ref_struct);
  SECTION( "Can identify number of matches" ) {
    REQUIRE(cohort.match_is_rare(0,A) == true);
    REQUIRE(cohort.match_is_rare(1,A) == false);
    REQUIRE(cohort.match_is_rare(2,A) == false);
  }
  // TODO deprecate
  SECTION( "Can extract the set of less-common alleles (as vector)" ) {
    vector<size_t> active_at_0 = cohort.get_active_rows(0, A);
    vector<size_t> active_at_1 = cohort.get_active_rows(1, A);
    vector<size_t> active_at_2 = cohort.get_active_rows(2, A);
    REQUIRE(active_at_0.size() == 1);
    REQUIRE(active_at_1.size() == 0);
    REQUIRE(active_at_2.size() == 2);
    REQUIRE(active_at_0[0] == 0);
    REQUIRE(active_at_2[0] == 2);
    REQUIRE(active_at_2[1] == 3);
  }
  SECTION( "Can extract the set of less-common alleles (as rowSet)" ) {
    rowSet active_at_0 = cohort.get_active_rowSet(0, A);
    rowSet active_at_1 = cohort.get_active_rowSet(1, A);
    rowSet active_at_2 = cohort.get_active_rowSet(2, A);
    REQUIRE(!active_at_0.empty());
    REQUIRE(active_at_1.empty());
    REQUIRE(*(active_at_0.begin()) == 0);
    rowSet::const_iterator it = active_at_2.begin();
    REQUIRE(*it == 2);
    ++it;
    REQUIRE(*it == 3);
  }
}

TEST_CASE( "haplotypeCohort construction", "[cohort][cohort-constructors]") {
  string ref_seq = "GATTACA";
  vector<size_t> positions = {1,4,5};
  siteIndex ref_struct = build_ref(ref_seq, positions);
  SECTION( "build-from-vectors" ) {
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
    
    ofstream testout;
    testout.open("testout.cohort", ios::out | ios::trunc);
    direct_cohort.serialize_human(testout);
    testout.close();
    REQUIRE(direct_cohort.get_n_sites() == 3);
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
  SECTION( "build from file" ) {
    ifstream testin;
    testin.open("testout.cohort", ios::in);
    siteIndex read_ref_struct(testin);
    haplotypeCohort read_cohort(testin, &read_ref_struct);
    testin.close();
    remove("testout.cohort");
    REQUIRE(read_cohort.get_n_sites() == 3);
    REQUIRE(read_cohort.number_matching(0,gap) == 2);
    REQUIRE(read_cohort.number_not_matching(0,gap) == 4);
  }   
  SECTION( "build-from-strings" ) {
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
    
    REQUIRE(string_cohort.get_n_sites() == 3);
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

TEST_CASE( "inputHaplotype", "[input-haplotype]" ) {
  // indices        0123456
  // sites           x  xx
  string ref_seq = "GATTACA";
  vector<size_t> positions = {1,4,5};
  siteIndex ref_struct = build_ref(ref_seq, positions);
  
  SECTION( "novel_SNVs are indexed as desired " ) {
    vector<alleleValue> alleles = {A, A};
    vector<size_t> counts = {1,2,3};
    inputHaplotype query = inputHaplotype(alleles, counts);
    REQUIRE(query.get_n_novel_SNVs(-1) == 1);
    REQUIRE(query.get_n_novel_SNVs(1) == 3);
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
    inputHaplotype query_0_string = inputHaplotype(seq_1.c_str(), ref_seq.c_str(), &ref_struct);
    inputHaplotype query_1_string = inputHaplotype(seq_2.c_str(), ref_seq.c_str(), &ref_struct);
    
    REQUIRE(query_0_direct.get_n_novel_SNVs(-1) ==
              query_0_string.get_n_novel_SNVs(-1));
    REQUIRE(query_0_direct.get_n_novel_SNVs(0) == 
              query_0_string.get_n_novel_SNVs(0));
    REQUIRE(query_1_direct.get_n_novel_SNVs(-1) == 
              query_1_string.get_n_novel_SNVs(-1));
    REQUIRE(query_1_direct.get_n_novel_SNVs(0) == 
              query_1_string.get_n_novel_SNVs(0));
    REQUIRE(query_0_direct.get_allele(0) == query_0_string.get_allele(0));
    REQUIRE(query_1_direct.get_allele(0) == query_1_string.get_allele(0));
  }
}

TEST_CASE( "penaltySet math", "[math]" ) {
  SECTION( "log-sum-exp implementation" ) {
    vector<double> R = {-1, -2, -3, -4};
    double naive_sum = logsum(logsum(-1, -2), logsum(-3, -4));
    double log_sum_exp = log_big_sum(R);
    REQUIRE(naive_sum == Approx(log_sum_exp));
    // repeated values
    R = {-1, -1, -2, -2};
    naive_sum = logsum(logsum(-1, -1), logsum(-2, -2));
    log_sum_exp = log_big_sum(R);
    REQUIRE(naive_sum == Approx(log_sum_exp));
    // similar numbers
    R = {-1.0000000001, -1, -2.000000001, -2};
    naive_sum = logsum(logsum(-1.000000001, -1), logsum(-2.000000001, -2));
    log_sum_exp = log_big_sum(R);
    REQUIRE(naive_sum == Approx(log_sum_exp));
  }
}

TEST_CASE( "Haplotype probabilities initial site", "[probability][probability-initial]" ) {
  SECTION( "Partial likelihoods at an initial site" ) {
    penaltySet penalties = penaltySet(-6, -9, 2);

    string ref_seq = "A";
    vector<size_t> positions = {0};
    siteIndex ref_struct = build_ref(ref_seq, positions);

    vector<string> allA = vector<string>(2, "A");
    vector<string> allT = vector<string>(2, "T");
    vector<string> AT = {"A","T"};    

    haplotypeCohort cohort_allA = haplotypeCohort(allA, &ref_struct);
    haplotypeCohort cohort_allT = haplotypeCohort(allT, &ref_struct);
    haplotypeCohort cohort_AT = haplotypeCohort(AT, &ref_struct);

    inputHaplotype query = inputHaplotype("A", ref_seq.c_str(), &ref_struct);

    fastFwdAlgState matrix_allA = fastFwdAlgState(&ref_struct, &penalties, 
                &cohort_allA);
    double probability_allA = matrix_allA.calculate_probability(&query);

    fastFwdAlgState matrix_allT = fastFwdAlgState(&ref_struct, &penalties, 
                &cohort_allT);
    double probability_allT = matrix_allT.calculate_probability(&query);

    fastFwdAlgState matrix_AT = fastFwdAlgState(&ref_struct, &penalties, 
                &cohort_AT);
    double probability_AT = matrix_AT.calculate_probability(&query);

    double expected_R0_in_allA_cohort = log(0.5) + penalties.one_minus_mu;
    double expected_R0_in_allT_cohort = log(0.5) + penalties.mu;
    double expected_R0_in_mixed_cohort = expected_R0_in_allA_cohort;
    double expected_R1_in_mixed_cohort = expected_R0_in_allT_cohort;
    double expected_probability_in_allA_cohort = log(2) 
              + expected_R0_in_allA_cohort;
    double expected_probability_in_allT_cohort = log(2) 
              + expected_R0_in_allT_cohort;
    double expected_probability_in_mixed_cohort = 
              logsum(expected_R0_in_mixed_cohort,expected_R1_in_mixed_cohort);

    matrix_allA.take_snapshot();
    matrix_allT.take_snapshot();
    matrix_AT.take_snapshot();

    REQUIRE(matrix_allA.R[0] == Approx(expected_R0_in_allA_cohort));
    REQUIRE(matrix_allT.R[0] == Approx(expected_R0_in_allT_cohort));

    REQUIRE(matrix_AT.R[0] == Approx(expected_R0_in_mixed_cohort));
    REQUIRE(matrix_AT.R[1] == Approx(expected_R1_in_mixed_cohort));

    REQUIRE(probability_allA == Approx(expected_probability_in_allA_cohort));
    REQUIRE(probability_allT == Approx(expected_probability_in_allT_cohort));
    REQUIRE(probability_AT == Approx(expected_probability_in_mixed_cohort));
  }
}

TEST_CASE( "Haplotype probabilities second site", "[probability][probability-site]" ) {
  SECTION( "Partial likelihoods at a second site" ) {
    penaltySet penalties = penaltySet(-6, -9, 4);
    
    string ref_seq = "AA";
    vector<size_t> positions = {0,1};
    siteIndex ref_struct = build_ref(ref_seq, positions);
    
    vector<string> allA = {"AA","TA","AA","TA"};
    vector<string> allT = {"AT","TT","AT","TT"};
    vector<string> AT = {"AA","TA","AT","TT"};    
    
    haplotypeCohort cohort_allA = haplotypeCohort(allA, &ref_struct);
    haplotypeCohort cohort_allT = haplotypeCohort(allT, &ref_struct);
    haplotypeCohort cohort_AT = haplotypeCohort(AT, &ref_struct);
    
    inputHaplotype query = inputHaplotype("AA", ref_seq.c_str(), &ref_struct);
    // fastFwdAlgState matrix_allA = fastFwdAlgState(&ref_struct, &penalties, 
    //             &cohort_allA);
    // fastFwdAlgState matrix_allT = fastFwdAlgState(&ref_struct, &penalties, 
    //             &cohort_allT);
    // fastFwdAlgState matrix_AT = fastFwdAlgState(&ref_struct, &penalties, 
    //             &cohort_AT);
    double R0A = log(0.25) + penalties.one_minus_mu;
    double R0T = log(0.25) + penalties.mu;
    double S0 = log(0.5) + logsum(penalties.mu, penalties.one_minus_mu);
    double ftR0A = penalties.R_coefficient + R0A;
    double ftR0T = penalties.R_coefficient + R0T;
    double pS = penalties.rho + S0;
    
    double expected_R0_in_allA_cohort = penalties.one_minus_mu 
              + logsum(ftR0A, pS);
    double expected_R1_in_allA_cohort = penalties.one_minus_mu 
              + logsum(ftR0T, pS);
    double expected_R0_in_allT_cohort = penalties.mu + logsum(ftR0A, pS);
    double expected_R1_in_allT_cohort = penalties.mu + logsum(ftR0T, pS);
    double expected_R0_in_mixed_cohort = expected_R0_in_allA_cohort;
    double expected_R1_in_mixed_cohort = expected_R1_in_allA_cohort;
    double expected_R2_in_mixed_cohort = expected_R0_in_allT_cohort;
    double expected_R3_in_mixed_cohort = expected_R1_in_allT_cohort;
    double expected_probability_in_allA_cohort = log(2) 
              + logsum(expected_R0_in_allA_cohort, expected_R1_in_allA_cohort);
    double expected_probability_in_allT_cohort = log(2) 
              + logsum(expected_R0_in_allT_cohort, expected_R1_in_allT_cohort);
    double expected_probability_in_mixed_cohort = 
              logsum(logsum(expected_R0_in_mixed_cohort, 
                            expected_R1_in_mixed_cohort), 
                     logsum(expected_R2_in_mixed_cohort, 
                            expected_R3_in_mixed_cohort));
      
    // matrix_allA.initialize_probability(&query);
    // matrix_allA.extend_probability_at_site(&query, 1);
    // matrix_allA.take_snapshot();
    // REQUIRE(matrix_allA.R[0] == Approx(expected_R0_in_allA_cohort));
    // REQUIRE(matrix_allA.R[1] == Approx(expected_R1_in_allA_cohort));
    // 
    // matrix_allT.initialize_probability(&query);
    // matrix_allT.extend_probability_at_site(&query, 1);
    // matrix_allT.take_snapshot();
    // REQUIRE(matrix_allT.R[0] == Approx(expected_R0_in_allT_cohort));
    // REQUIRE(matrix_allT.R[1] == Approx(expected_R1_in_allT_cohort));
    // 
    // matrix_AT.initialize_probability(&query);
    // matrix_AT.extend_probability_at_site(&query, 1);
    // matrix_AT.take_snapshot();
    // REQUIRE(matrix_AT.R[0] == Approx(expected_R0_in_mixed_cohort));
    // REQUIRE(matrix_AT.R[1] == Approx(expected_R1_in_mixed_cohort));
    // REQUIRE(matrix_AT.R[2] == Approx(expected_R2_in_mixed_cohort));
    // REQUIRE(matrix_AT.R[3] == Approx(expected_R3_in_mixed_cohort));
    // 
    // REQUIRE(matrix_allA.S == Approx(expected_probability_in_allA_cohort));
    // REQUIRE(matrix_allT.S == Approx(expected_probability_in_allT_cohort));
    // REQUIRE(matrix_AT.S == Approx(expected_probability_in_mixed_cohort));
  }
}

TEST_CASE( "Haplotype probabilities at terminal span", "[probability][probability-terminal-span]" ) {
  SECTION( "Partial likelihoods at a span following a site" ) {
    penaltySet penalties = penaltySet(-6, -9, 3);
    
    string ref_seq = "AAAAAA";
    vector<size_t> positions = {0};
    siteIndex ref_struct = build_ref(ref_seq, positions);
    
    vector<string> haplotypes = {
      "AAAAAA",
      "AAAAAA",
      "TAAAAA"
    };
    haplotypeCohort cohort = haplotypeCohort(haplotypes, &ref_struct);
    
    inputHaplotype query_0_aug = inputHaplotype("AAAAAA", ref_seq.c_str(),
              &ref_struct);
    inputHaplotype query_1_aug = inputHaplotype("AAAAAT", ref_seq.c_str(), 
              &ref_struct);
    
    double mu = penalties.mu;
    double mu_c = penalties.one_minus_mu;
    double alpha = penalties.R_coefficient;
    double rho = penalties.rho;
    double alpha_l = 5 * alpha;
    
    double lR_i = mu_c - log(3);
    double lR_im = mu - log(3);
    double lS_i = logsum(lR_i + log(2), lR_im);
    
    double RHS = lS_i - log(3) + log1p(-exp(alpha_l));
    
    double expected_R0_for_0_mismatch_haplotype = 
              mu_c*5 + logsum(alpha_l + lR_i, RHS);
    double expected_R2_for_0_mismatch_haplotype = 
              mu_c*5 + logsum(alpha_l + lR_im, RHS);
    double expected_R0_for_1_mismatch_haplotype =
              mu_c*4 + mu + logsum(alpha_l + lR_i, RHS);
    double expected_R2_for_1_mismatch_haplotype =
              mu_c*4 + mu + logsum(alpha_l + lR_im, RHS);
    double expected_probability_for_0_mismatch_haplotype =
              mu_c*5 + lS_i;
    double expected_probability_for_1_mismatch_haplotype =
              mu_c*4 + mu + lS_i;
    
    fastFwdAlgState matrix_0_aug = fastFwdAlgState(&ref_struct, &penalties, 
                &cohort);
    matrix_0_aug.initialize_probability(&query_0_aug);
    matrix_0_aug.extend_probability_at_span_after(&query_0_aug, 0);
    matrix_0_aug.take_snapshot();

    REQUIRE(matrix_0_aug.S == Approx(expected_probability_for_0_mismatch_haplotype));
    REQUIRE(matrix_0_aug.R[0] == Approx(expected_R0_for_0_mismatch_haplotype));
    REQUIRE(matrix_0_aug.R[2] == Approx(expected_R2_for_0_mismatch_haplotype));

    fastFwdAlgState matrix_1_aug = fastFwdAlgState(&ref_struct, &penalties, 
                &cohort);
    matrix_1_aug.initialize_probability(&query_1_aug);
    matrix_1_aug.extend_probability_at_span_after(&query_1_aug, 0);
    matrix_1_aug.take_snapshot();
  
    REQUIRE(matrix_1_aug.R[0] == Approx(expected_R0_for_1_mismatch_haplotype));
    REQUIRE(matrix_1_aug.R[2] == Approx(expected_R2_for_1_mismatch_haplotype));  
    REQUIRE(matrix_1_aug.S == Approx(expected_probability_for_1_mismatch_haplotype));
  }
}
  
TEST_CASE( "Haplotype probabilities at initial span", "[probability][probability-initial-span]" ) {
  SECTION( "Partial likelihoods at an initial span" ) {
    penaltySet penalties = penaltySet(-6, -9, 3);
    
    string ref_seq = "AAAAAA";
    vector<size_t> positions = {5};
    siteIndex ref_struct = build_ref(ref_seq, positions);
    
    vector<string> haplotypes = {
      "AAAAAA",
      "AAAAAA",
      "AAAAAT"
    };
    haplotypeCohort cohort = haplotypeCohort(haplotypes, &ref_struct);
    
    inputHaplotype query_0_aug = 
              inputHaplotype("AAAAAA", ref_seq.c_str(), &ref_struct);
    inputHaplotype query_1_aug = 
              inputHaplotype("TAAAAA", ref_seq.c_str(), &ref_struct);
    
    fastFwdAlgState matrix_0_aug = fastFwdAlgState(&ref_struct, &penalties, 
                &cohort);
    fastFwdAlgState matrix_1_aug = fastFwdAlgState(&ref_struct, &penalties, 
                &cohort);
    
    double mu = penalties.mu;
    double mu_c = penalties.one_minus_mu;
    
    double expected_R0_for_0_mismatch_haplotype = 6*mu_c - log(3);
    double expected_R2_for_0_mismatch_haplotype = 5*mu_c + mu - log(3);
    double expected_R0_for_1_mismatch_haplotype = 5*mu_c + mu - log(3);
    double expected_R2_for_1_mismatch_haplotype = 4*mu_c + 2*mu - log(3);
    double expected_probability_for_0_mismatch_haplotype = 5*mu_c + logdiff(log(2), log(7) + mu) - log(3);
    double expected_probability_for_1_mismatch_haplotype = 4*mu_c + mu + logdiff(log(2), log(7) + mu) - log(3);
    
    matrix_0_aug.initialize_probability(&query_0_aug);
    matrix_0_aug.take_snapshot();
    REQUIRE(matrix_0_aug.R[0] == Approx(expected_R0_for_0_mismatch_haplotype));
    REQUIRE(matrix_0_aug.R[2] == Approx(expected_R2_for_0_mismatch_haplotype));
    REQUIRE(matrix_0_aug.S == Approx(expected_probability_for_0_mismatch_haplotype));
    
    matrix_1_aug.initialize_probability(&query_1_aug);
    matrix_1_aug.take_snapshot();
    REQUIRE(matrix_1_aug.R[0] == Approx(expected_R0_for_1_mismatch_haplotype));
    REQUIRE(matrix_1_aug.R[2] == Approx(expected_R2_for_1_mismatch_haplotype));
    REQUIRE(matrix_1_aug.S == Approx(expected_probability_for_1_mismatch_haplotype));
  }
}

// TEST_CASE( "Haplotype probabilities at series of sites", "[probability][probability-multiple-sites]" ) {
//   SECTION( "Partial likelihoods at a series of sites" ) {
//     penaltySet penalties = penaltySet(-6, -9, 3);
//     double mu_c = penalties.one_minus_mu;
//     double mu = penalties.mu;
//     double rho = penalties.rho;
//     double alpha = penalties.R_coefficient;
//     double l2_mu = logdiff(log(2), log(7) + mu);
//     
//     string ref_seq = "AAAAA";
//     vector<string> haplotypes = {
//       "AAAAA",
//       "TATAT",
//       "ATACA"
//     };
//     string query_string = "ATAAT";
//     
//     // Matrix of cohort-haplotype matching to query. 1 indicates match; 0
//     // indicates non-match
//     // 10110
//     // 00011
//     // 11100
//     
//     vector<size_t> positions = {0,1,2,3,4};
//     siteIndex ref_struct = build_ref(ref_seq, positions);
//     haplotypeCohort cohort = haplotypeCohort(haplotypes, &ref_struct);
//     inputHaplotype query = inputHaplotype(query_string.c_str(), ref_seq.c_str(), &ref_struct);
//     fastFwdAlgState matrix = fastFwdAlgState(&ref_struct, &penalties, 
//                 &cohort);
//       
//     vector<double> Ss;
//     vector<vector<double> > Rs;
//     
//     // R0 -- initial site [o|x|o]
//     Rs.push_back({
//         -log(3) + mu_c,
//         -log(3) + mu,
//         -log(3) + mu_c
//     });
//     // S0
//     Ss.push_back(-log(3) + l2_mu);
//     
//     // R1 [x|x|o]
//     Rs.push_back({
//         mu + logsum(alpha + Rs[0][0], rho + Ss[0]),
//         mu + logsum(alpha + Rs[0][1], rho + Ss[0]),
//         mu_c + logsum(alpha + Rs[0][2], rho + Ss[0])
//     });
//     // S1
//     Ss.push_back(logsum(Rs[1][0], logsum(Rs[1][1],Rs[1][2])));
//     
//     // R2 [o|x|o]
//     Rs.push_back({
//         mu_c + logsum(alpha + Rs[1][0], rho + Ss[1]),
//         mu + logsum(alpha + Rs[1][1], rho + Ss[1]),
//         mu_c + logsum(alpha + Rs[1][2], rho + Ss[1])
//     });
//     // S2
//     Ss.push_back(logsum(Rs[2][0], logsum(Rs[2][1],Rs[2][2])));
//     
//     // R3 [o|o|x]
//     Rs.push_back({
//         mu_c + logsum(alpha + Rs[2][0], rho + Ss[2]),
//         mu_c + logsum(alpha + Rs[2][1], rho + Ss[2]),
//         mu + logsum(alpha + Rs[2][2], rho + Ss[2])
//     });
//     // S3
//     Ss.push_back(logsum(Rs[3][0], logsum(Rs[3][1],Rs[3][2])));
//     
//     // R4 [x|o|x]
//     Rs.push_back({
//         mu + logsum(alpha + Rs[3][0], rho + Ss[3]),
//         mu_c + logsum(alpha + Rs[3][1], rho + Ss[3]),
//         mu + logsum(alpha + Rs[3][2], rho + Ss[3])
//     });
//     // S4
//     Ss.push_back(logsum(Rs[4][0], logsum(Rs[4][1],Rs[4][2])));
//     
//     // 0 match / 1 not / 2 match
//     // all active since initial site
//     matrix.initialize_probability(&query);
//     REQUIRE(matrix.R[0] == Approx(Rs[0][0]));
//     REQUIRE(matrix.R[1] == Approx(Rs[0][1])); 
//     REQUIRE(matrix.R[2] == Approx(Rs[0][2])); 
//     REQUIRE(matrix.S == Approx(Ss[0]));
//     
//     DPUpdateMap current_map;
//     bool match_is_rare;
//     
//     // 0 not / 1 not / 2 match
//     // non-match common / active = {2}
//     match_is_rare = cohort.match_is_rare(1, T);
//     REQUIRE(match_is_rare == true);
// 
//     current_map = penalties.get_current_map(matrix.S, match_is_rare);
//     REQUIRE(current_map.coefficient == Approx((mu + alpha)));
//     REQUIRE(current_map.constant == Approx((rho + Ss[0] - alpha)));
//     
//     matrix.extend_probability_at_site(&query, 1);
//     REQUIRE(matrix.R[2] == Approx(Rs[1][2])); 
//     REQUIRE(matrix.S == Approx(Ss[1]));
//     
//     // 0 match / 1 not / 2 match
//     // match common / active = {1}
//     match_is_rare = cohort.match_is_rare(2, A);
//     REQUIRE(match_is_rare == false);
// 
//     current_map = penalties.get_current_map(matrix.S, match_is_rare);
//     REQUIRE(current_map.coefficient == Approx((mu_c + alpha)));
//     REQUIRE(current_map.constant == Approx((rho + Ss[1] - alpha)));
//         
//     matrix.extend_probability_at_site(&query, 2);
//     REQUIRE(matrix.R[1] == Approx(Rs[2][1])); 
//     REQUIRE(matrix.S == Approx(Ss[2]));
// 
//     // 0 match / 1 match / 2 not
//     // match common / active = {2}
//     matrix.extend_probability_at_site(&query, 3);
//     REQUIRE(matrix.R[2] == Approx(Rs[3][2])); 
//     REQUIRE(matrix.S == Approx(Ss[3]));
// 
//     // 0 not / 1 match / 2 not
//     // not-match common / active = {1}
//     matrix.extend_probability_at_site(&query, 4);
//     matrix.take_snapshot();
//     REQUIRE(matrix.R[0] == Approx(Rs[4][0]));
//     REQUIRE(matrix.R[1] == Approx(Rs[4][1])); 
//     REQUIRE(matrix.R[2] == Approx(Rs[4][2]));
//     REQUIRE(matrix.S == Approx(Ss[4]));
//   }
// }

TEST_CASE( "Haplotype probabilities at homogeneous sites", "[probability][probability-homogeneous-sites]" ) {
  SECTION( "Partial likelihood at a series of sites without variation equals the likelihood at a span of equivalent length" ) {
    penaltySet penalties = penaltySet(-6, -9, 3);
    
    string ref_seq = "AAAAA";
    vector<size_t> positions_span = {0};
    siteIndex ref_struct_span = build_ref(ref_seq, 
                positions_span);

    vector<size_t> positions = {0,1,2,3,4};
    siteIndex ref_struct = build_ref(ref_seq, positions);
        
    vector<string> haplotypes = {
      "AAAAA",
      "AAAAA",
      "TAAAA"
    };

    haplotypeCohort cohort = haplotypeCohort(haplotypes, &ref_struct);
    haplotypeCohort cohort_span = haplotypeCohort(haplotypes, &ref_struct_span);    
    
    inputHaplotype query = inputHaplotype(ref_seq.c_str(), ref_seq.c_str(), &ref_struct);
    inputHaplotype query_span = inputHaplotype(ref_seq.c_str(), ref_seq.c_str(), 
                &ref_struct_span);
    
    fastFwdAlgState matrix = fastFwdAlgState(&ref_struct, &penalties, 
                &cohort);
    fastFwdAlgState matrix_span = fastFwdAlgState(&ref_struct_span, &penalties, 
                &cohort_span);
        
    double probability = matrix.calculate_probability(&query);
    double probability_span = matrix_span.calculate_probability(&query_span);

    REQUIRE(matrix.R[0] == Approx(matrix_span.R[0]));
    REQUIRE(matrix.R[1] == Approx(matrix_span.R[1])); 
    REQUIRE(matrix.R[2] == Approx(matrix_span.R[2])); 
    REQUIRE(probability == Approx(probability_span));
  }
}

TEST_CASE( "Fast method gives same result as classical", "[probability][sublinear-correctness]" ) {
  SECTION( "Fast, linear and quadratic forward algorithms give same result " ) {
    penaltySet penalties = penaltySet(-6, -9, 5);
    vector<size_t> positions = {0, 1, 2, 3, 4};
    vector<vector<alleleValue> > haplotypes = {
      {A, A, A, T, G},
      {A, C, A, A, G},
      {C, A, A, T, G},
      {A, C, C, A, G},
      {A, A, A, A, G},
    };
    siteIndex reference(positions, 5);
    haplotypeCohort cohort(haplotypes, &reference);
    
    slowFwdSolver linear_fwd(&reference, &penalties, &cohort);
    slowFwdSolver quadratic_fwd(&reference, &penalties, &cohort);
    fastFwdAlgState fast_fwd(&reference, &penalties, &cohort);
    
    vector<alleleValue> query = {A, A, A, A, A};
    inputHaplotype query_ih(query, vector<size_t>(5, 0), &reference, 0, 5);
    
    double result_fast = fast_fwd.calculate_probability(&query_ih);
    double result_linear = linear_fwd.calculate_probability_linear(query, 0);
    double result_quad = quadratic_fwd.calculate_probability_quadratic(query, 0);
    REQUIRE(result_fast == Approx(result_linear));
    REQUIRE(result_quad == Approx(result_linear));
  }
}

// TEST_CASE( "Relative indexing works", "[haplotype][reference][input]" ) {
//   //                01234567890123456789
//   // sites              4    9    4
//   string ref_seq = "AAAAAAAAAAAAAAAAAAAA";
//   // left            AAA
//   // middle              AAA 
//   // right                         AAA    
//   // 9 14                   AAAAAAAA
//   //                01234567890123456789
//   
//   vector<size_t> positions = {4, 9, 14};
//   siteIndex ref_struct = build_ref(ref_seq, positions);
//   
//   SECTION( "relative indices are produced as desired " ) {
//     string eight = "TTTTTTTT";
//     string three = "AAA";
//     string six = "TTTTTT";
//     inputHaplotype nosites_left = inputHaplotype(three.c_str(), ref_seq.c_str(), 
//                                                  &ref_struct, 1, 3);
//     inputHaplotype nosites_middle = inputHaplotype(three.c_str(), ref_seq.c_str(),
//                                                   &ref_struct, 5, 3);
//     inputHaplotype nosites_right = inputHaplotype(three.c_str(), ref_seq.c_str(),
//                                                     &ref_struct, 15, 3);
//     inputHaplotype overlap_9_16 = inputHaplotype(eight.c_str(), ref_seq.c_str(),
//                                                  &ref_struct, 8, 8);
//     inputHaplotype overlap_4_9 = inputHaplotype(six.c_str(), ref_seq.c_str(),
//                                                  &ref_struct, 4, 6);
//     REQUIRE(nosites_left.has_sites() == false);
//     REQUIRE(nosites_right.has_sites() == false);
//     REQUIRE(nosites_middle.has_sites() == false);
//     REQUIRE(overlap_9_16.has_sites() == true);
//     REQUIRE(overlap_9_16.number_of_sites() == 2);
//     REQUIRE(overlap_9_16.get_left_tail() == 1);
//     REQUIRE(overlap_9_16.get_site_index(0) == 1);
//     REQUIRE(overlap_9_16.get_n_novel_SNVs(-1) == 1);
//     REQUIRE(overlap_9_16.get_n_novel_SNVs(0) == 4);
//     REQUIRE(overlap_9_16.get_n_novel_SNVs(1) == 1);
//     REQUIRE(overlap_4_9.has_left_tail() == false);
//     REQUIRE(overlap_4_9.number_of_sites() == 2);
//     REQUIRE(overlap_4_9.has_span_after(1) == false);
//     REQUIRE(overlap_9_16.get_n_novel_SNVs(0) == 4);
//   }
// }

TEST_CASE( "Delay map object works ", "[delay][DPMap]" ) {
  SECTION( "DPMaps work" ) {
    DPUpdateMap ID = DPUpdateMap(0);
    DPUpdateMap scale = DPUpdateMap(-2);
    REQUIRE(scale.is_degenerate());
    DPUpdateMap test1 = DPUpdateMap(-3, -4);
    DPUpdateMap test2 = DPUpdateMap(-5, -6);
    REQUIRE(ID.of(ID) == ID);
    REQUIRE(ID.of(scale) == scale);
    REQUIRE(scale.of(ID) == scale);
    REQUIRE(ID.of(-1.0) == -1.0);
    REQUIRE(ID.of(test1) == test1);
    REQUIRE(test1.of(ID) == test1);
    REQUIRE(scale.of(scale) == DPUpdateMap(-4.0));
    REQUIRE(scale.of(-1.0) == Approx(-3.0));
    REQUIRE((test1.of(scale)).coefficient == Approx(-5.0));
    REQUIRE((test1.of(scale)).constant == Approx(-2.0));
    REQUIRE((scale.of(test1)).coefficient == Approx(-5.0));
    REQUIRE((scale.of(test1)).constant == Approx(-4.0));
    REQUIRE((test1.of(test2)).coefficient == Approx(-8.0));
    REQUIRE((test1.of(test2)).constant == Approx(logsum(-6, 1)));
    DPUpdateMap test3 = test1.compose(test2);
    double expected_coefficient = -3 - 5;
    double expected_constant = logsum(-6, -4 + 5);
    REQUIRE(test3.coefficient == Approx(expected_coefficient));
    REQUIRE(test3.constant == Approx(expected_constant));
  }
}

TEST_CASE( "Delay map structure is built correctly ", "[delay][lazyEval-build]" ) {
  SECTION( "Building and accessing delayedEvalMaps" ) {
    delayedEvalMap map = delayedEvalMap(5);
    REQUIRE(map.get_map_indices().size() == 5);
    REQUIRE(map.get_map(0).is_identity() == true);
    // Can we add and then read a value?
    // Add a map with value (1.0, 1.0)
    map.add_eqclass(DPUpdateMap(1.0, 1.0));
    // Assign row 0 to map 1.0 at eqclass-index 0
    map.remove_row_from_eqclass(0);
    map.assign_row_to_newest_eqclass(0);
    REQUIRE(map.get_coefficient(0) == Approx(1));
    // Can we add and then read a second value?
    // Add a map with value (2.0, 2.0)
    map.add_eqclass(DPUpdateMap(2.0, 2.0));
    // Assign row 1 to map 2.0 at eqclass-index 1
    map.remove_row_from_eqclass(1);
    map.assign_row_to_newest_eqclass(1);
    REQUIRE(map.get_coefficient(1) == Approx(2));
    // Can we overload a eqclass with two rows?
    // Assign row 2 to map 2.0 at eqclass-index 1
    map.remove_row_from_eqclass(2);
    map.assign_row_to_newest_eqclass(2);
    REQUIRE(map.get_coefficient(2) == Approx(2));
    // Can we remove a row from a eqclass?
    map.remove_row_from_eqclass(1);
    REQUIRE(map.get_map_indices()[1] == 5);
    // When we empty all contents of a eqclass, do we return it to the set of eqclasses
    // which can be filled?
    // Remove row 1 from its eqclass 0. This should also delete eqclass 1 since it is
    // now empty
    map.remove_row_from_eqclass(0);
    // This new map should go in eqclass-index 1, which was just emptied
    map.add_eqclass(DPUpdateMap(3.0, 3.0));
    // Add a new row to this map
    map.assign_row_to_newest_eqclass(3);
    // This row should get eqclass-index 0
    REQUIRE(map.get_map_indices()[3] == 1);
    REQUIRE(map.get_coefficient(3) == Approx(3.0));
  }
}

TEST_CASE( "Delay map structure stores values correctly ", "[delay][storage]" ) {
  SECTION( "Updating maps performs correct arithmetic" ) {
    delayedEvalMap map = delayedEvalMap(3);
    map.stage_map_for_site(DPUpdateMap(-2.0, -3.0));
    REQUIRE(map.get_map_history()[0].is_identity());
    REQUIRE(map.get_map_history()[1] == DPUpdateMap(-2.0, -3.0));
    REQUIRE(map.row_updated_to(0) == 0);
    map.hard_update_all();
    // map.update_maps({0});
    REQUIRE(map.row_updated_to(0) == 1);
    double m_of_one = -2.0 + logsum(-3.0, -1.0);
    REQUIRE(map.get_coefficient(0) == -2.0);
    REQUIRE(map.get_constant(0) == -3.0);
    REQUIRE(map.get_map(0).of(-1.0) == m_of_one);
  }
}

TEST_CASE( "Delay map structure stores history correctly ", "[delay][encode-history]" ) {
  SECTION( "We can track sites and times-of-update in delayedEvalMaps" ) {
    delayedEvalMap map2 = delayedEvalMap(2);
    map2.add_eqclass(DPUpdateMap(-1.0, -1.0));
    map2.assign_row_to_newest_eqclass(0);
    REQUIRE(map2.last_update(0) == 0);
    map2.stage_map_for_site(DPUpdateMap(-3.0, -3.0));
    map2.hard_update_all();
    REQUIRE(map2.last_update(0) == 1);
    REQUIRE(map2.last_update(1) == 1);
  }
}

TEST_CASE( "PenaltySet gives right values", "[penaltyset]") {
  double eps = 0.0000001;
  penaltySet penalties = penaltySet(-6, -9, 3);
  double S = -1;
  double expected;
  double one_minus_mu = penalties.one_minus_mu;
  double mu = penalties.mu;
  double rho = penalties.rho;
  double alpha = penalties.R_coefficient;
  double one_minus_2mu = penalties.one_minus_2mu;
  
  expected = mu + alpha;
  REQUIRE(penalties.get_current_map(S, true).coefficient == Approx(expected));
  expected = rho - alpha + S;
  REQUIRE(penalties.get_current_map(S, true).constant == Approx(expected));
  expected = one_minus_mu + alpha;
  REQUIRE(penalties.get_current_map(S, false).coefficient == Approx(expected));
  expected = rho - alpha + S;
  REQUIRE(penalties.get_current_map(S, false).constant == Approx(expected));
  
  vector<double> summands = {-2, -3};
  expected = logsum(mu + S, one_minus_2mu - one_minus_mu + logsum(-2, -3));
  penalties.update_S(S, summands, true);
  REQUIRE(S == Approx(expected));
  expected = logdiff(one_minus_mu + S, one_minus_2mu - mu + logsum(-2, -3));
  penalties.update_S(S, summands, false);
  REQUIRE(S == Approx(expected));
}

// TEST_CASE( "Delay-maps perform correct state-update calculations", "[delay][delay-probability]" ) {
//   double eps = 0.0000001;
//   penaltySet penalties = penaltySet(-6, -9, 3);
//   
//   string ref_seq = "AAAAA";
//   vector<size_t> positions = {0,1,2,3,4};
//   siteIndex ref_struct = build_ref(ref_seq, positions);
//   
//   vector<string> haplotypes = {
//     "AAAAA",
//     "TATAT",
//     "ATACA"
//   };
//   // 10110
//   // 00011
//   // 11100
//   haplotypeCohort cohort = haplotypeCohort(haplotypes, &ref_struct);
//   
//   string query_string = "ATAAT";
//   inputHaplotype query = inputHaplotype(query_string.c_str(), ref_seq.c_str(), &ref_struct);
//   
//   fastFwdAlgState matrix = fastFwdAlgState(&ref_struct, &penalties, 
//               &cohort);
//   
//   double mu_c = penalties.one_minus_mu;
//   double mu = penalties.mu;
//   double rho = penalties.rho;
//   double alpha = penalties.R_coefficient;
//   double S;
//   double mu2 = penalties.one_minus_2mu;
  // SECTION( "Delay map maintains correct states for all rows in matrix" ) {
  //   fastFwdAlgState matrix = fastFwdAlgState(&ref_struct, &penalties, 
  //               &cohort);
  //   matrix.initialize_probability(&query);
  //   REQUIRE(matrix.get_maps().number_of_eqclasses() == 1);
  //   REQUIRE(matrix.get_maps().row_updated_to(0) == 0);
  //   REQUIRE(matrix.get_maps().get_map(0).is_identity() == true);
  //   REQUIRE(matrix.get_maps().get_map(1).is_identity() == true);
  //   REQUIRE(matrix.get_maps().get_map(2).is_identity() == true);
  //   matrix.extend_probability_at_site(&query, 1);
  //   REQUIRE(matrix.get_maps().number_of_eqclasses() == 2);
  //   REQUIRE(matrix.get_maps().row_updated_to(0) == 1);
  //   REQUIRE(matrix.get_maps().row_updated_to(2) == 1);
  //   REQUIRE(matrix.get_maps().get_map(0).is_identity() == false);
  //   REQUIRE(matrix.get_maps().get_map(1).is_identity() == false);
  //   REQUIRE(matrix.get_maps().get_map(2).is_identity() == true);
  //   matrix.extend_probability_at_site(&query, 2);
  //   REQUIRE(matrix.get_maps().number_of_eqclasses() == 3);
  //   REQUIRE(matrix.get_maps().row_updated_to(0) == 2);
  //   REQUIRE(matrix.get_maps().row_updated_to(1) == 2);
  //   REQUIRE(matrix.get_maps().row_updated_to(2) == 1);
  //   REQUIRE(matrix.get_maps().get_map(0).is_identity() == false);
  //   REQUIRE(matrix.get_maps().get_map(1).is_identity() == true);
  //   REQUIRE(matrix.get_maps().get_map(2).is_identity() == true);
  //   matrix.extend_probability_at_site(&query, 3);
  //   REQUIRE(matrix.get_maps().number_of_eqclasses() == 3);
  //   REQUIRE(matrix.get_maps().row_updated_to(0) == 2);
  //   REQUIRE(matrix.get_maps().row_updated_to(1) == 2);
  //   REQUIRE(matrix.get_maps().row_updated_to(2) == 3);
  //   REQUIRE(matrix.get_maps().get_map(0).is_identity() == false);
  //   REQUIRE(matrix.get_maps().get_map(1).is_identity() == true);
  //   REQUIRE(matrix.get_maps().get_map(2).is_identity() == true);
  //   matrix.extend_probability_at_site(&query, 4);
  //   REQUIRE(matrix.get_maps().number_of_eqclasses() == 3);
  //   REQUIRE(matrix.get_maps().row_updated_to(0) == 2);
  //   REQUIRE(matrix.get_maps().row_updated_to(1) == 4);
  //   REQUIRE(matrix.get_maps().row_updated_to(2) == 3);    
  //   REQUIRE(matrix.get_maps().get_map(0).is_identity() == false);
  //   REQUIRE(matrix.get_maps().get_map(1).is_identity() == true);
  //   REQUIRE(matrix.get_maps().get_map(2).is_identity() == true);
  //   matrix.get_maps().hard_update_all();
  //   REQUIRE(matrix.get_maps().row_updated_to(0) == 4);
  //   REQUIRE(matrix.get_maps().row_updated_to(1) == 4);
  //   REQUIRE(matrix.get_maps().row_updated_to(2) == 4);
  //   REQUIRE(matrix.get_maps().get_map(0).is_identity() == false);
  //   REQUIRE(matrix.get_maps().get_map(1).is_identity() == true);
  //   REQUIRE(matrix.get_maps().get_map(2).is_identity() == false);
  //   matrix.get_maps().hard_clear_all();
  //   REQUIRE(matrix.get_maps().number_of_eqclasses() == 1);
  //   REQUIRE(matrix.get_maps().get_map(0).is_identity() == true);
  //   REQUIRE(matrix.get_maps().get_map(1).is_identity() == true);
  //   REQUIRE(matrix.get_maps().get_map(2).is_identity() == true);
  // }
//   SECTION( "Partial likelihood is correctly calculated at a series of sites" ) {
//     vector<double> correct_coefficients;
//     vector<vector<double> > correct_constants;
//     // 10110
//     // 00011
//     // 11100
//     // @ 0, haplotypes 0 and 2 match, match more common
//     // work done at all since initial
//     matrix.initialize_probability(&query);
//     REQUIRE(matrix.get_maps().get_map(0).is_identity());
//     REQUIRE(matrix.get_maps().get_map(1).is_identity());
//     REQUIRE(matrix.get_maps().get_map(2).is_identity());
//     // work done at 2; match is rare
//     matrix.extend_probability_at_site(&query, 1);
//     S = logdiff(log(2), log(7) + mu) - log(3);
//     DPUpdateMap I = DPUpdateMap(0);
//     DPUpdateMap m1 = penalties.get_non_match_map(S);
//     REQUIRE(matrix.get_maps().get_map(0).constant == Approx(m1.constant));
//     REQUIRE(matrix.get_maps().get_map(0).coefficient == Approx(m1.coefficient));
//     REQUIRE(matrix.get_maps().get_map(1).constant == Approx(m1.constant));
//     REQUIRE(matrix.get_maps().get_map(1).coefficient == Approx(m1.coefficient));
//     REQUIRE(matrix.get_maps().get_map(2).is_identity());
//     // work done at 1; non-match is rare
//     matrix.extend_probability_at_site(&query, 2);
//     S = logsum(mu + S, mu2 + logsum(alpha + mu_c - log(3), rho + S));
//     DPUpdateMap M2 = penalties.get_match_map(S);
//     REQUIRE(matrix.get_maps().get_map(0).constant == Approx((M2.of(m1)).constant));
//     REQUIRE(matrix.get_maps().get_map(0).coefficient == Approx((M2.of(m1)).coefficient));
//     REQUIRE(matrix.get_maps().get_map(1).is_identity());
//     REQUIRE(matrix.get_maps().get_map(2).is_identity());
//     S = matrix.S;
//     matrix.extend_probability_at_site(&query, 3);
//     DPUpdateMap M3 = penalties.get_match_map(S);
//     REQUIRE(matrix.get_maps().get_map(0).constant == Approx((M2.of(m1)).constant));
//     REQUIRE(matrix.get_maps().get_map(0).coefficient == Approx((M2.of(m1)).coefficient));
//     REQUIRE(matrix.get_maps().get_map(1).is_identity());
//     REQUIRE(matrix.get_maps().get_map(2).is_identity());
//     S = matrix.S;
//     matrix.extend_probability_at_site(&query, 4);
//     DPUpdateMap m4 = penalties.get_non_match_map(S);
//     REQUIRE(matrix.get_maps().get_map(0).constant == Approx((M2.of(m1)).constant));
//     REQUIRE(matrix.get_maps().get_map(0).coefficient == Approx((M2.of(m1)).coefficient));
//     REQUIRE(matrix.get_maps().get_map(1).is_identity());
//     REQUIRE(matrix.get_maps().get_map(2).is_identity());
//     matrix.take_snapshot();
//     REQUIRE(matrix.get_maps().get_map(0).is_identity());
//     REQUIRE(matrix.get_maps().get_map(1).is_identity());
//     REQUIRE(matrix.get_maps().get_map(2).is_identity());
//   }
// }