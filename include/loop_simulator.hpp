#ifndef INCLUDE_LOOP_SIMULATOR_HPP
#define INCLUDE_LOOP_SIMULATOR_HPP

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <deque>
#include <random>
#include <algorithm>
#include <loop.hpp>
#include <bond_array.hpp>

class loop_simulator
{
public:

  // constructor and destructor
  loop_simulator();
  ~loop_simulator();

  void prng_seed(int s);

  void initialize_loop_simulator(double basal_death_prob,
        double step_prob,
        double stall_death_prob,
        double knockoff,
        double bypass,
        int N, // size of chromosome
        int numSmc,
        int smcWidth);

  int wrapped_index(int pos);
  void occupy_around(int pos);
  void release_around(int pos);

  void birth(int i);
  void death();
  void step(bool tag_extruded);
  void steps(int N_steps, bool tag_extruded);

  std::vector<int> neighbors(int pos);
  bool is_empty_around(int pos);
  bool no_forks_around(int pos);
  int get_next_position(int pos, int direction);
  int get_next_vacancy(int pos, int direction);

  // read and write the state
  void write_state(std::string st_filename);
  void read_state(std::string st_filename);

  void set_N(int N_new);
  void set_M(int M_new);

  // WCM patch: numSmc is authoritative from loop_params (not replication-scaled).
  int get_numSmc_initial() const { return numSmc_initial; }
  int get_loaded_loop_count() const { return static_cast<int>(loops1.size()); }
  void sync_M_to_loaded();

  int read_loop_params(std::string loop_param_filename);

  std::vector<bond> get_loop_bonds();
  std::vector<loop> get_loops();
  std::vector<std::deque<int>> get_extruded();

private:

  // loop topology
  // std::vector<loop> loops;
  std::mt19937 rand_eng;

  int N; // Current length of the chromosome (increases during replication)
  int N_initial; // Initial chromosome length (before replication starts)
  int fork_width;
  int left_fork;
  int right_fork;
  int M; // Number of SMCs
  int numSmc_initial; // WCM patch: from loop_params numSmc (not replication-scaled)

  double falloff;
  double step_prob;
  double stall_falloff;
  double knockoff;
  double bypass;

  std::vector<int> loops1;
  std::vector<int> loops2;
  std::vector<int> stalled1;
  std::vector<int> stalled2;
  std::vector<int> occupied;

  std::vector<std::deque<int>> extruded;

  int smc_width;

};

#endif
