#include <iostream>
#include <fstream>
#include <string>
#include <sys/time.h>
#include <cmath>
#include "probability.hpp"
#include "reference.hpp"

// PROFILER: this application is used for profiling the time complexity of our sublinear Li Stephens
// application

namespace profiler_tools{
  bool validate_penalties(double mutation, double recombination) {
    return (mutation > 0 && recombination > 0);
  }
  
  // bool validate_sizes(const char* max, const char* number) {
  //   
  // }
  
  vector<size_t> get_sizes(bool logarithmic, size_t max, size_t number) {
    vector<size_t> to_return = vector<size_t>(number);
    if(!logarithmic) {
      double step_size = max / number;
      for(size_t i = 0; i < number; i++) {
        to_return[i] = (size_t)((i + 1) * step_size);
      }
    } else {
      double step_size = log(max) / number;
      for(size_t i = 0; i < number; i++) {
        to_return[i] = (size_t)(exp((i + 1) * step_size));
      }
    }
    return to_return;
  }
  
  // bool validate_lengths(const char* min, const char* max, const char* number) {
  //   
  // }
  
  vector<size_t> get_lengths(bool logarithmic, size_t min, size_t max, size_t number) {
    vector<size_t> to_return = vector<size_t>(number);
    if(!logarithmic) {
      double step_size = (max - min)/(number - 1);
      for(size_t i = 0; i < number; i++) {
        to_return[i] = min + (size_t)(i * step_size);
      }
    } else {
      double step_size = (log(max) - log(min)) / (number - 1);
      for(size_t i = 0; i < number; i++) {
        to_return[i] = (size_t)(min * exp(i * step_size));
      }
    }
    return to_return;
  }
  
  bool size_compatible(const haplotypeCohort* cohort, const vector<size_t>& sizes) {
    size_t number_of_haplotypes = cohort->get_n_haplotypes();
    for(size_t i = 0; i < number_of_haplotypes; i++) {
      if(sizes[i] == 0 || sizes[i] > number_of_haplotypes) {
        return false;
      }
    }
    return true;
  }
  
  bool length_compatible(const siteIndex* reference, const vector<size_t>& lengths) {
    size_t number_of_sites = reference->number_of_sites();
    for(size_t i = 0; i < number_of_sites; i++) {
      if(lengths[i] == 0 || lengths[i] > number_of_sites) {
        return false;
      }
    }
    return true;
  }
  
  void write_experiment(ostream& out, const vector<size_t>& sizes, const vector<size_t>& lengths, size_t replicates, double mutation_pen, double recombination_pen) {
    out << sizes.size() << "\t";
    out << lengths.size() << "\t";
    out << replicates << endl;
    for(size_t i = 0; i < sizes.size(); i++) {
      out << sizes[i] << "\t";
    }
    for(size_t i = 0; i < lengths.size(); i++) {
      out << lengths[i] << "\t";
    }
    out << endl << mutation_pen << "\t" << recombination_pen;
  }
  
  void read_experiment(istream& in, vector<size_t>& sizes, vector<size_t>& lengths, size_t& replicates, double& mutation_pen, double& recombination_pen) {
    size_t n_sizes;
    size_t n_lengths;
    in >> n_sizes;
    in >> n_lengths;
    in >> replicates;
    sizes.resize(n_sizes);
    lengths.resize(n_lengths);
    for(size_t i = 0; i < n_sizes; i++) {
      in >> sizes[i];
    }
    for(size_t i = 0; i < n_lengths; i++) {
      in >> lengths[i];
    }
    in >> mutation_pen;
    in >> recombination_pen;
  }
  
  
  inputHaplotype* random_haplo(haplotypeCohort* cohort, siteIndex* reference, size_t generations, penaltySet* penalties) {
    size_t start_site = 0;
    if(start_site < cohort->get_n_sites()) {
      size_t end_site = cohort->get_n_sites() - 1;
      vector<alleleValue> random_haplo = cohort->rand_desc_haplo(generations, penalties->rho, penalties->mu, start_site, end_site);
      inputHaplotype* to_return = new inputHaplotype(random_haplo, vector<size_t>(random_haplo.size(), 0), reference, cohort->get_reference()->get_position(start_site), cohort->get_reference()->get_position(end_site) - cohort->get_reference()->get_position(start_site));
      to_return->validate();
      return to_return;
    } else {
      return new inputHaplotype(reference);
    }
  }
}

int main(int argc, char* argv[]) {
	// -- input handling -------------------------------------------------------------------------------------------------
  bool print_commands = false;
  bool write_experiment = false;
  
  string out_path;
  string slls_path;
  string experiment_path;
  bool print_cout = false;
  
  bool log_spacing;
  vector<size_t> lengths;
  vector<size_t> sizes;
  size_t replicates;  

  double mutation_penalty;
  double recombination_penalty;
  
  size_t haplotype_generator_generations = 3;
  size_t replicates_per_datapoint = 3;
  size_t LINEAR_MAX_SAMPLES = 10000;
  size_t QUADRATIC_MAX_NK2 = 10000000;
  
  if(argc == 1 || (argc == 2 && (argv[1] == string("help") || argv[1] == string("commands")))) {
    print_commands = true;
  } else if((argc == 13 || argc == 14 || argc == 17 || argc == 18) && argv[1] == string("define-experiment")) {
    write_experiment = true;
    
    size_t command_pos = 2;
    bool sizes_built = false;
    bool lengths_built = false;
    bool replicates_built = false;
    bool recombinations_built = false;
    bool mutations_built = false;
    bool print_define_experiment_commands = false;
    
    if(argc == 13 || argc == 17) {
      print_cout = true;
    } else {
      command_pos = 3;
      experiment_path = argv[2];
    }
    if(argc == 13 || argc == 14) {
      mutation_penalty = 2.3 * 9;
      recombination_penalty = 2.3 * 6;
    }
    
    while(!(sizes_built && lengths_built && replicates_built) || ((argc < 17) || !(recombinations_built && mutations_built))) {
      bool repeated = false;
      if(argv[command_pos] == string("sizes")) {
        if(sizes_built) {
          repeated = true;
        } else {
          if(argv[command_pos + 1] == string("log") || argv[command_pos + 1] == string("logarithmic")) {
            log_spacing = true;
          } else if(argv[command_pos + 1] == string("lin") || argv[command_pos + 1] == string("linear")) {
            log_spacing = false;
          } else {
            print_define_experiment_commands = true;
            break;
          }
          size_t par_max_size = atol(argv[command_pos + 2]);
          size_t par_size_steps = atol(argv[command_pos + 3]);
          command_pos += 4;
          sizes = profiler_tools::get_sizes(log_spacing, par_max_size, par_size_steps);
          sizes_built = true;
        }
      } else if(argv[command_pos] == string("lengths")) {
        if(lengths_built) {
          repeated = true;
        } else {
          if(argv[command_pos + 1] == string("log") || argv[command_pos + 1] == string("logarithmic")) {
            log_spacing = true;
          } else if(argv[command_pos + 1] == string("lin") || argv[command_pos + 1] == string("linear")) {
            log_spacing = false;
          } else {
            print_define_experiment_commands = true;
            break;
          }
          size_t par_min_length = atol(argv[command_pos + 2]);
          size_t par_max_length = atol(argv[command_pos + 3]);
          size_t par_length_steps = atol(argv[command_pos + 4]);
          command_pos += 5;
          lengths = profiler_tools::get_lengths(log_spacing, par_min_length, par_max_length, par_length_steps);
          lengths_built = true;
        }
      } else if(argv[command_pos] == string("replicates")) {
        if(replicates_built) {
          repeated = true;
        } else {
          replicates = atol(argv[command_pos + 1]);
          command_pos += 2;
          replicates_built = true;
        }
      } else if(argc == 17 || argc == 18) {
        if(argv[command_pos] == string("recombination")) {
          if(recombinations_built) {
            repeated = true;
          } else {
            recombination_penalty = atof(argv[command_pos + 1]);
            command_pos += 2;
            recombinations_built = true;
          }
        } else if(argv[command_pos] == string("mutation")) {
          if(mutations_built) {
            repeated = true;
          } else {
            mutation_penalty = atof(argv[command_pos + 1]);
            command_pos += 2;
            mutations_built = true;
          }
        } else {
          print_define_experiment_commands = true;
          break;
        }
      } else {
        print_define_experiment_commands = true;
        break;
      }
      if(repeated) {
        cerr << "parameter " << argv[command_pos] << " repeated" << endl;
        return 1;
      }
    }
    
    if(print_define_experiment_commands) {
      cerr << "-- command options are ------------------------------------" << endl;
      cerr << "profiler define-experiment [experiment_path] sizes [log OR lin] [max] [number] lengths [log OR lin] [min] [max] [number] replicates [number]" << endl;
      cerr << "\t with additional optional parameters" << endl;
      cerr << "\t mutation [log-penalty] recombination [log-penalty]" << endl;
      cerr << "\t (defaults mutation 9, recombination 6)" << endl;
      return 1;
    }
  } else if(argv[1] == string("speed")) {
    if(argc == 4 || argc == 5) {
      slls_path = argv[2];
      experiment_path = argv[3];
      if(argc == 4) {
        print_cout = true;
      } else {
        out_path = argv[4];
      }
    } else {
      cerr << "-- command options are ------------------------------------" << endl;
      cerr << "profiler speed [slls_path] [experiment_path] [output-path]" << endl;
      cerr << "profiler deltas [slls_path] [output-path]" << endl;
      cerr << "profiler space [slls_path] [output-path]" << endl;
      cerr << "\t if output-path is not defined, defaults to cout" << endl;
      return 1;
    }
  } else if(argv[1] == string("deltas")) {
    if(argc == 3) {
      print_cout = true;
    } else if(argc == 4) {
      out_path = argv[3];
    } else {
      cerr << "-- command options are ------------------------------------" << endl;
      cerr << "profiler deltas [slls_path] [output-path]" << endl;
      cerr << "\t if output-path is not defined, defaults to cout" << endl;
      return 1;
    }
  } else if(argv[1] == string("space")) {
    if(argc == 3) {
      print_cout = true;
    } else if(argc == 4) {
      out_path = argv[3];
    } else {
      cerr << "-- command options are ------------------------------------" << endl;
      cerr << "profiler space [slls_path] [output-path]" << endl;
      cerr << "\t if output-path is not defined, defaults to cout" << endl;
      return 1;
    }
  } else {
    print_commands = true;
  }  
  
  if(print_commands) {
    cerr << "-- command options are --------------------------------" << endl;
    cerr << "define-experiment : writes an experiment parameter file" << endl;
    cerr << "speed \t : runs speed profiling" << endl;
    cerr << "deltas \t : runs delta-parameter profiling" << endl;
    cerr << "space \t : runs space profiling" << endl << endl;
    cerr << "-- where parameters are -------------------------------" << endl;
    cerr << "profiler define-experiment [experiment_path] sizes [log OR lin] [max] [number] lengths [log OR lin] [min] [max] [number] replicates [number]" << endl;
    cerr << "\t with additional optional parameters" << endl;
    cerr << "\t mutation [log-penalty] recombination [log-penalty]" << endl;
    cerr << "\t (defaults mutation 9, recombination 6)" << endl;
    cerr << "profiler speed [slls_path] [experiment_path] [output-path]" << endl;
    cerr << "profiler deltas [slls_path] [output-path]" << endl;
    cerr << "profiler space [slls_path] [output-path]" << endl;
    cerr << "\t if output-path is not defined, defaults to cout" << endl;
    return 1;
  }
  
  if(write_experiment) {
    ofstream expt_pars_out;
    if(!print_cout) {
      expt_pars_out.open(experiment_path, ios::out | ios::trunc);
      profiler_tools::write_experiment(expt_pars_out, sizes, lengths, replicates, recombination_penalty, mutation_penalty);
      expt_pars_out.close();
    } else {
      profiler_tools::write_experiment(cout, sizes, lengths, replicates, recombination_penalty, mutation_penalty);
    }
    return 0;
  }
  
  ifstream slls_in;
  slls_in.open(slls_path, ios::in);
  if(!slls_in.is_open()) {
    cerr << "failed to open slls index" << endl;
    return 1;
  }  
  cerr << "loading reference from slls index " << slls_path << endl;
  siteIndex* reference = new siteIndex(slls_in);
  cerr << "loading haplotype cohort from slls index" << endl;
  haplotypeCohort* cohort = new haplotypeCohort(slls_in, reference);
  slls_in.close();
  size_t cohort_size = cohort->get_n_haplotypes();
  size_t n_trials = atoi(argv[2]);
  
  ofstream file_output;
  file_output.open(out_path, ios::out | ios::trunc);
  ostream& output = print_cout ? cout : file_output;
  
  // build penalty container
  penaltySet* penalties = new penaltySet(recombination_penalty, mutation_penalty, cohort_size);
    
  for(size_t j = 0; j < n_trials; j++) {
    for(size_t n_bp_it = 0; n_bp_it < lengths.size(); n_bp_it++) {
      for(size_t k_it = 0; k_it < sizes.size(); ) {
        size_t k = sizes[k_it];
        size_t n_bp = lengths[n_bp_it];
        size_t interval_start = reference->rand_interval_start(n_bp);
        size_t interval_end = interval_start + n_bp;
        size_t site_start = reference->find_site_above(interval_start);
        size_t site_end = reference->find_site_below(interval_end);
        haplotypeCohort* new_cohort = cohort->subset(site_start, site_end, k);
        siteIndex* new_reference = cohort->get_reference();
        
        if(new_reference->number_of_sites() > 1) {
        	inputHaplotype* query_ih = profiler_tools::random_haplo(new_cohort, new_reference, haplotype_generator_generations, penalties);
          
          if(query_ih->number_of_sites() > 1) {
            cerr << "iteration "<< k_it << " of " << sizes.size() << " size, " << n_bp_it << " of " << lengths.size() << " length, " << j << " of " << replicates << " replicates " << endl;
            ++k_it;
            size_t start_site = query_ih->get_start_site();
            size_t end_site = start_site + query_ih->number_of_sites() - 1;    
            
            double time_used_fast = 0;
          	double time_used_quad = 0;
          	double time_used_linear = 0;
                      	
          	for(size_t i = 0; i < replicates_per_datapoint; i++) {
              fastFwdAlgState* haplotype_matrix = new fastFwdAlgState(new_reference, penalties, new_cohort);
            	slowFwdSolver* linear_fwd = new slowFwdSolver(new_reference, penalties, new_cohort);
            	slowFwdSolver* quadratic_fwd = new slowFwdSolver(new_reference, penalties, new_cohort);

            	struct timeval tv1, tv2, tv3, tv4;
            	gettimeofday(&tv1, NULL);
            	double result = haplotype_matrix->calculate_probability(query_ih);
              gettimeofday(&tv2, NULL);
              time_used_fast += (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 + (double) (tv2.tv_sec - tv1.tv_sec);

              double result_slowq;
              if(k * k * new_reference->number_of_sites() > QUADRATIC_MAX_NK2) {
                time_used_quad = 0;
              } else {
                gettimeofday(&tv2, NULL);
                result_slowq = quadratic_fwd->calculate_probability_quadratic(query_ih);
                gettimeofday(&tv3, NULL);
          	    time_used_quad += (double) (tv3.tv_usec - tv2.tv_usec) / 1000000 + (double) (tv3.tv_sec - tv2.tv_sec);
              }
              
              double result_slowl;
              if(cohort_size > LINEAR_MAX_SAMPLES) {
                time_used_linear = 0;
              } else {
                gettimeofday(&tv3, NULL);
              	result_slowl = linear_fwd->calculate_probability_linear(query_ih);
                gettimeofday(&tv4, NULL);
                time_used_linear += (double) (tv4.tv_usec - tv3.tv_usec) / 1000000 + (double) (tv4.tv_sec - tv3.tv_sec);
              }
              
              cerr << result << " ?= " << result_slowq << " ?= " << result_slowl << endl;
              
              delete haplotype_matrix;
            	delete quadratic_fwd;
            	delete linear_fwd;
            }
            
            size_t sum_of_information_content = new_cohort->sum_information_content(query_ih->get_alleles(), query_ih->get_start_site());
            
            output << "time fast\t" << time_used_fast/replicates_per_datapoint << "\ttime linear\t" << time_used_linear/replicates_per_datapoint << "\ttime quadratic\t" << time_used_quad/replicates_per_datapoint << "\tsum of information_content\t" << sum_of_information_content << "\tsites in haplotype\t" << end_site - start_site << "\tcohort size\t" << k << "\tregion length in bp\t" << n_bp << endl;
          }
          delete query_ih;
        }
        delete new_cohort;
        delete new_reference;
      }
    }
  }
  delete cohort;
  delete reference;
  delete penalties;
  if(print_cout) {
    file_output.close();
  }
  return 0;
}