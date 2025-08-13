#include <loop_simulator.hpp>

// constructor
loop_simulator::loop_simulator()
{
  rand_eng.seed(0);
  // 100 beads in 1 s
  // 1 bead in 10 ms
  // basal death prob: 0.0004 per second, or 2500 seconds lifetime
  // step prob: 1
  // stall death prob: 0.001 to 0.005 per second
  // knockoff prob: ? set to 1-bypass
  // bypass prob: lifetime of ~20 s: 0.03 to 0.05 per second
  // initialize_loop_simulator(0.000004,1.0,0.00005,0.9995,0.0005,54338,50,1);
  initialize_loop_simulator(0.00005,1.0,0.00005,0.0,1.0,54338,50,1);
}


// destructor
loop_simulator::~loop_simulator()
{
}


void loop_simulator::prng_seed(int s)
{
    rand_eng.seed(s);
}

void loop_simulator::initialize_loop_simulator(
        double basal_death_prob,
        double step_prob,
        double stall_death_prob,
        double knockoff,
        double bypass,
        int N, // size of chromosome
        int numSmc,
        int smcWidth
    )
{
        this->N = N;
        int midpoint = 27169; // 54338 is even so technically 27169 and 27170 both straddle origin (to start)
        this->fork_width = (N - 54338) / 2;
        this->left_fork = midpoint - fork_width;
        this->right_fork = midpoint + fork_width + 1;

        this->M = numSmc;

        this->falloff = basal_death_prob;
        this->step_prob = step_prob;
        this->stall_falloff = stall_death_prob;
        this->knockoff = knockoff;
        this->bypass = bypass; // rate of bypassing

        this->loops1 = std::vector<int>(this->M, 0);
        this->loops2 = std::vector<int>(this->M, 0);
        this->occupied = std::vector<int>(this->N + 1, 0);
        this->stalled1 = std::vector<int>(this->M, 0);
        this->stalled2 = std::vector<int>(this->M, 0);

        this->smc_width = smcWidth; // if smc_width=1, we occupy_around pos-1, pos and pos+1

        for (int ind = 0; ind < this->M; ++ind) {
          this->birth(ind);
        }

        this->extruded = std::vector<std::deque<int>>(this->N);
}

int loop_simulator::wrapped_index(int idx) {
    if (idx < 1) return N + idx;
    if (idx > N) return idx - N;
    return idx;
}

void loop_simulator::occupy_around(int pos) {
    for (int j = -smc_width; j <= smc_width; ++j) {
        occupied[wrapped_index(pos + j)]++; // need2b careful of wraparound logic?
    }
}

void loop_simulator::release_around(int pos) {
    for (int j = -smc_width; j <= smc_width; ++j) {
        occupied[wrapped_index(pos + j)]--; // need2b careful of wraparound logic?
    }
}

std::vector<int> loop_simulator::neighbors(int pos) {
    std::vector<int> result;
    for (int j = -smc_width; j <= smc_width; ++j) {
        result.push_back(occupied[wrapped_index(pos + j)]); // need2b careful of wraparound logic?
    }
    return result;
}

bool loop_simulator::is_empty_around(int pos) {
    for (int i : neighbors(pos)) {
        if (i > 0) {
            return false;
        }
    }
    return true;
}

bool loop_simulator::no_forks_around(int pos) {
    if (this->fork_width == 0) return true ;

    for (int j = -smc_width; j <= smc_width; ++j) {
        if (wrapped_index(pos + j) == left_fork || wrapped_index(pos + j) == right_fork) {
            return false;
        }
    }
    return true;
}

void loop_simulator::birth(int i) {

    while (true) {
        std::uniform_int_distribution<int> unif_dist(1, N);
        int pos = unif_dist(rand_eng);

        if (pos <= 54338) {
            // Chromosome 1 (original / mother strand or mother-indexed daughter)

            int idx1 = pos;
            int idx2 = ((pos + 1 + 2*smc_width - 1) % 54338) + 1;  // wrap from 54338 back to 1

            std::cout << "SMC: idx1 = " << idx1 << ", idx2 = " << idx2 << ", distance = " << (idx2 > idx1 ? idx2 - idx1 : 54338 - idx1 + idx2) << "\n" << std::endl;

            // Reject if either spot is already occupied or if too close to replication forks
            if (!(is_empty_around(idx1)
                && is_empty_around(idx2)
                && no_forks_around(idx1)
                && no_forks_around(idx2))) continue;

            loops1[i] = idx1;
            loops2[i] = idx2;
            occupy_around(idx1);
            occupy_around(idx2);
            return;
        } else {
            // Chromosome 2 (replicated daughter strand)

            int idx1 = pos;
            int idx2 = pos + 1 + 2*smc_width;


            if (!(is_empty_around(idx1)
                && is_empty_around(idx2))
                && idx1 >= 54339 + smc_width
                && idx2 <= N - smc_width) continue;

            loops1[i] = idx1;
            loops2[i] = idx2;
            occupy_around(idx1);
            occupy_around(idx2);
            return;
        }
    }
}


void loop_simulator::death()
{
std::uniform_real_distribution<double> dist(0.0, 1.0);
        for (int i = 0; i < M; ++i) {

            // mark for dissociation if either side is stalled
            //if (stalled1[i] == 1 || stalled2[i] == 1) {
            //    stalled1[i] = 1;
            //    stalled2[i] = 1;
            //}

            double falloff1 = stalled1[i] == 0 ? falloff : stall_falloff;
            double falloff2 = stalled2[i] == 0 ? falloff : stall_falloff;
            double falloff_rate = std::max(falloff1, falloff2);
            if (dist(rand_eng) < falloff_rate) {
                release_around(loops1[i]);
                release_around(loops2[i]);
                stalled1[i] = 0;
                stalled2[i] = 0;
                birth(i);
            }
        }
}

int loop_simulator::get_next_position(int pos, int direction) {
    // get next position based on topology
   int new_pos;


   // for now skip trying to debug the above code
    if (pos <= 54338) {
        int region_size = 54338;
        int zero_based_pos = pos - 1;
        int wrapped = (zero_based_pos + direction) % region_size;
        if (wrapped < 0) wrapped += region_size;
        new_pos = wrapped + 1;
    } else {
        int region_start = 54339;
        int region_size = N - 54338;
        int zero_based_pos = pos - region_start;
        int wrapped = (zero_based_pos + direction) % region_size;
        if (wrapped < 0) wrapped += region_size;
        new_pos = wrapped + region_start;
    }

    return new_pos;
}

int loop_simulator::get_next_vacancy(int pos, int direction) {
    int candidate = pos;
    for (int i = 0; i<15; ++i) { // look 15 away
        candidate = get_next_position(candidate, direction);
        if (is_empty_around(candidate) && no_forks_around(candidate)) {
            return candidate;
        }
    }
    return pos; // we searched
}

void loop_simulator::step(bool tag_extruded)
{
        std::uniform_real_distribution<> dis(0.0, 1.0);

        for (int i = 0; i < M; ++i) {

            // update positions
            int cur1 = loops1[i];
            int cur2 = loops2[i];

            // if SMC position is in the pool of cytosolic SMCs, skip it
            if (cur1 >= 2 * N || cur2 >= 2 * N) continue;

            // new positions with the periodic boundaries
            int newCur1 = get_next_position(cur1, -1);
            int newCur2 = get_next_position(cur2, 1);

            //assert(newCur1 <= 2 * N - 1);
            //assert(newCur2 <= 2 * N - 1);
	        stalled1[i] = 0;
            stalled2[i] = 0;

            if (cur1 == newCur1) stalled1[i] = 1;  //shouldn't happen
            if (cur2 == newCur2) stalled2[i] = 1;

            //std::cout << "Moving SMC..." <<std::endl;
            //  move each side only if the other side is "free" or stall/knock off
            if (stalled1[i] == 0) {
		        //std::cout << "SMC is not stalled..." <<std::endl;
                if (dis(rand_eng) <= step_prob) {
		            //std::cout << "Attempting to take a step to newCur1 = " << newCur1  <<std::endl;
                    if (neighbors(newCur1).front() == 0 && no_forks_around(newCur1)) {
                        // take a normal step
			            // std::cout << "Taking a normal step..." <<std::endl;
                        occupy_around(newCur1);
                        release_around(cur1);
                        loops1[i] = newCur1;
                        if (tag_extruded) extruded[i].push_front(newCur1);
                    } else {
			            // std::cout << "Can't take a normal step, either stalling or bypassing..." <<std::endl;
                        double rate_sum = knockoff + bypass;
                        if (dis(rand_eng) <= rate_sum) {
                            if (dis(rand_eng) <= knockoff / rate_sum || occupied[newCur2] > 2) {
                                // stall and maybe dissociate
				                // std::cout << "Marking as stalled..." <<std::endl;
                                stalled1[i] = 1;
                            } else {
                                // One could imagine that bypassing is harder the more SMCs are in the traffic jam
                                // We can allow SMCs to bypass the entire jam with the requirement that they traverse
                                // no greater than 50 nm during the bypass
				                std::cout << "Bypassing..." <<std::endl;
				                // STILL NEED TO ADD LOGIC FOR BYPASSING REPLICATION FORK
				                newCur1 = get_next_vacancy(cur1, -1);
                                occupy_around(newCur1);
                                release_around(cur1);
                                loops1[i] = newCur1;
                                if (cur1==newCur1) {
                                    std::cout << "... but there was no vacancy" <<std::endl;
                                } else {
                                    if (tag_extruded) extruded[i].push_front(newCur1);
                                }
                            }
                        }
                    }
                }
            }

            if (stalled2[i] == 0) {
                if (dis(rand_eng) <= step_prob) {
                    if (neighbors(newCur2).back() == 0 && no_forks_around(newCur2)) {
                        // take a normal step
                        occupy_around(newCur2);
                        release_around(cur2);
                        loops2[i] = newCur2;
                        if (tag_extruded) extruded[i].push_back(newCur2);
                    } else {
                        double rate_sum = knockoff + bypass;
                        if (dis(rand_eng) <= rate_sum) {
                            if (dis(rand_eng) <= knockoff / rate_sum || occupied[newCur2] > 2) {
                                // stall and maybe dissociate
                                stalled2[i] = 1;
                            } else {
                                // bypass
                                newCur2 = get_next_vacancy(cur2, 1);
                                occupy_around(newCur2);
                                release_around(cur2);
                                loops2[i] = newCur2;
                                if (cur2==newCur2) {
                                    std::cout << "... but there was no vacancy" <<std::endl;
                                } else {
                                    if (tag_extruded) extruded[i].push_back(newCur2);
                                }
                            }
                        }
                    }
                }
            }


        }
}

void loop_simulator::steps(int N_steps, bool tag_extruded)
{
    // reset extruded
    this->extruded = std::vector<std::deque<int>>(this->N);

    for (int i = 0; i < N_steps; ++i) {
            if (i% 1000 == 0) std::cout << i << std::endl;
            death();
            step(tag_extruded);
    }
}

void loop_simulator::write_state(std::string st_filename) {

    std::fstream st_file;

    // Open file in append mode
    st_file.open(st_filename, std::ios::out | std::ios::app);

    if (!st_file.is_open()) {
        std::cout << "ERROR: file not opened in write_state" << std::endl;
    } else {
        // Write a separator header indicating the number of loops
        st_file << "Number of loops: " << loops1.size() << "\n";

        int midpoint = 27169;
        int fork_width = (N - 54338) / 2;
        int left_fork = midpoint - fork_width;
        int right_fork = midpoint + fork_width + 1;
        st_file << "Replication forks: " << left_fork << ", " << right_fork << "\n";

        for (size_t i = 0; i < loops1.size(); ++i) {
            st_file << loops1[i] << "\t" << loops2[i] << "\t"
                    << stalled1[i] << "\t" << stalled2[i];

            // Write neighborhood around loops1[i]
            st_file << "\t[";
            for (int offset = -2*smc_width; offset <= 2*smc_width; ++offset) {
                int idx = wrapped_index(loops1[i] + offset);
                st_file << occupied[idx];
                if (offset < 2*smc_width) st_file << ",";
            }
            st_file << "]";

            // Write neighborhood around loops2[i]
            st_file << "\t[";
            for (int offset = -2*smc_width; offset <= 2*smc_width; ++offset) {
                int idx = wrapped_index(loops2[i] + offset);
                st_file << occupied[idx];
                if (offset < 2*smc_width) st_file << ",";
            }
            st_file << "]";

            st_file << "\n";
        }

        // Add a newline for separation between different runs
        st_file << "\n";
    }

    st_file.close();
}

void loop_simulator::read_state(std::string st_filename) {
    std::fstream st_file;
    st_file.open(st_filename, std::ios::in);
    std::cout << "Clearing old loops..." << std::endl;
    loops1.clear();
    loops2.clear();

    std::string line;
    std::cout << "Reading in new loops..." << std::endl;
    while (std::getline(st_file, line)) {
        // Skip metadata lines
        if (line.find("Number of loops") != std::string::npos ||
            line.find("Replication forks") != std::string::npos ||
            line.empty())
            continue;

        std::istringstream iss(line);
        int val1, val2;
        if (iss >> val1 >> val2) {
            loops1.push_back(val1);
            loops2.push_back(val2);
            occupy_around(val1);
            occupy_around(val2);
        }
    }

    st_file.close();
}

std::vector<bond> loop_simulator::get_loop_bonds() {
    bond loop_bond;
    std::vector<bond> loop_bonds;
    int id = N + 1;

    std::vector<loop> loops = get_loops();

    for (loop l : loops)
    {
        if (l.get_a_bound())
        {
            loop_bond.id = id;
            loop_bond.type = 2;
            loop_bond.i = l.get_a();
            loop_bond.j = l.get_h();
            loop_bonds.push_back(loop_bond);
            id += 1;
        }
    }
    return loop_bonds;
}

std::vector<loop> loop_simulator::get_loops() {
    std::vector<loop> loops;

    for (size_t i = 0; i < loops1.size(); ++i) {
        loop loopi(0,0);
        loopi.set_a(loops1[i]);
        loopi.set_h(loops2[i]);
        loopi.set_a_bound(true);
        loopi.set_h_bound(true);
        loops.push_back(loopi);
    }

    return loops;
}

std::vector<std::deque<int>> loop_simulator::get_extruded() {
    return extruded;
}

void loop_simulator::set_N(int N_new) {
    int N_old = this->N;  // store the old value
    std::vector<int> old_occupied = this->occupied;  // copy old occupied array

    this->N = N_new;
    std::cout << "Resizing occupied array" << std::endl;
    this->occupied = std::vector<int>(this->N + 1, 0);  // resize to new length

    std::cout << "Copying elements from old occupied" << std::endl;
    int min_len = std::min(N_old + 1, N_new + 1);  // Avoid out-of-bounds
    for (int i = 0; i < min_len; ++i) {
        this->occupied[i] = old_occupied[i];
    }

    std::cout << "Inferring fork positions from new chromosome length" << std::endl;
    this->fork_width = (this->N - 54338) / 2;
    int midpoint = 27169; // 54338 is even so technically 27169 and 27170 both straddle origin (to start)
    this->left_fork = midpoint - fork_width;
    this->right_fork = midpoint + fork_width + 1;

    std::cout << "Done setting N" << std::endl;

    std::cout << "Re-indexing SMC on right daughter, and knocking off all SMC that are too close to or straddle a replisome..." << std::endl;
    for (int i = 0; i < M; ++i) {
        if ((loops1[i] >= left_fork && loops1[i] <= right_fork && loops2[i] >= right_fork) //    L 1 R 2
            || (loops1[i] <= left_fork && loops2[i] >= left_fork  && loops2[i] <= right_fork) // 1 L 2 R
            || (loops1[i] >= left_fork && loops1[i] <= right_fork && loops2[i] <= left_fork)  // 2 L 1 R
            || (loops1[i] >= right_fork && loops2[i] >= left_fork  && loops2[i] <= right_fork) // L 2 R 1
            || !no_forks_around(loops1[i])
            || !no_forks_around(loops2[i])
        ) {
            std::cout << "Knocked off SMC " << i << "! Binding elsewhere..." << std::endl;
            release_around(loops1[i]);
            release_around(loops2[i]);
            stalled1[i] = 0;
            stalled2[i] = 0;
            birth(i);
        }
        if (loops1[i] > 54338 && loops2[i] > 54338) { // 
            release_around(loops1[i]);
            release_around(loops2[i]);
            stalled1[i] = 0;
            stalled2[i] = 0;
            loops1[i] = loops1[i] + (N_new - N_old) / 2;
            loops2[i] = loops2[i] + (N_new - N_old) / 2;
            occupy_around(loops1[i]);
            occupy_around(loops2[i]);
        }
    }


}

void loop_simulator::set_M(int M_new) {
    int M_old = this->M;
    this->M = M_new;

    // Resize vectors while preserving old values
    std::cout << "Resizing loop arrays" << std::endl;
    this->loops1.resize(this->M, 0);
    this->loops2.resize(this->M, 0);
    this->stalled1.resize(this->M, 0);
    this->stalled2.resize(this->M, 0);
    std::cout << "Birthing new loops" << std::endl;

    // Birth only the new loops
    for (int ind = M_old; ind < this->M; ++ind) {
        this->birth(ind);
    }
}
